#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#if HIGH_PREC == GMP_HIGH_PREC
 #include <gmpxx.h>
#endif

#include "base/random.hpp"
#include "communicate/communicate.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp"
#include "hmc/multipseudo/multipseudo_rhmc_step.hpp"
#include "hmc/hmc.hpp"
#include "geometry/geometry_eo.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/rat_approx.hpp"
#include "new_types/su3.hpp"
#include "operations/remez/remez_algorithm.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //fourth root of 2, used to extend the range of eigenvalues
  const double enl_gen=pow(2,0.25);
  
  //Return the maximal eigenvalue of the staggered Dirac operator for the passed quark
  THREADABLE_FUNCTION_6ARG(max_eigenval, double*,eig_max, quark_content_t*,quark, quad_su3**,eo_conf, clover_term_t**,Cl, quad_u1**,backfield, int,niters)
  {
    pseudofermion_t in(quark->discretiz);
    pseudofermion_t temp1(quark->discretiz);
    pseudofermion_t temp2(quark->discretiz); //not used for stag...
    pseudofermion_t out(quark->discretiz);
    
    //generate the random field
    in.fill();
    
    //perform initial normalization
    double init_norm=in.normalize();
    verbosity_lv3_master_printf("Init norm: %lg\n",init_norm);
    
    //prepare the ingredients
    add_backfield_with_stagphases_to_conf(eo_conf,backfield);
    inv_clover_term_t *invCl_evn=NULL;
    if(ferm_discretiz::include_clover(quark->discretiz))
      {
	invCl_evn=nissa_malloc("invCl_evn",loc_volh,inv_clover_term_t);
	invert_twisted_clover_term(invCl_evn,quark->mass,quark->kappa,Cl[EVN]);
      }
    
    //apply the vector niter times normalizing at each iter
    int iter=0;
    int is_increasing=1;
    double old_eig_max;
    
    do
      {
	switch(quark->discretiz)
	  {
	  case ferm_discretiz::ROOT_STAG:
	    apply_stD2ee_m2(out.stag,eo_conf,temp1.stag,sqr(quark->mass),in.stag);break;
	  case ferm_discretiz::ROOT_TM_CLOV:
	    tmclovDkern_eoprec_square_eos(out.Wils,temp1.Wils,temp2.Wils,eo_conf,quark->kappa,Cl[ODD],invCl_evn,quark->mass,in.Wils);
	    break;
	  default:
	    crash("not supported yet");
	  }
	
	//compute the norm
	old_eig_max=*eig_max;
	(*eig_max)=in.normalize(out);
	
	if((iter++)>0) is_increasing=(*eig_max/old_eig_max-1>1e-14);
	verbosity_lv2_master_printf("max_eigen search mass %lg, iter %d, eig %16.16lg\n",quark->mass,iter,*eig_max);
      }
    while(iter<niters&&is_increasing);
    
    //remove the background field
    rem_backfield_with_stagphases_from_conf(eo_conf,backfield);
    
    //assume a 10% excess
    (*eig_max)*=1.1;
    verbosity_lv2_master_printf("max_eigen mass %lg: %16.16lg\n",quark->mass,*eig_max);
  }
  THREADABLE_FUNCTION_END
  
  //check that an approximation is valid in the interval passed
  bool check_approx_validity(rat_approx_t &appr,double eig_min,double eig_max,double expo,double maxerr)
  {
    //compute condition numbers
    double cond_numb=eig_max/eig_min;
    double cond_numb_stored=appr.maximum/appr.minimum;
    
    //check that approximation is valid
    bool valid=true;
    
    //check that has been allocated
    if(valid)
      {
	bool allocated=(appr.degree()!=0);
	if(!allocated) verbosity_lv2_master_printf(" Not allocated\n");
	valid&=allocated;
      }
    
    //check that exponent is the same
    if(valid)
      {
	double expo_stored=(double)appr.num/appr.den;
	bool expo_is_the_same=(fabs(expo-expo_stored)<1e-14);
	if(!expo_is_the_same) verbosity_lv2_master_printf(" Approx stored is for x^%lg, required is for x^%lg\n",expo_stored,expo);
	valid&=expo_is_the_same;
      }
    
    //check that it can be fitted
    if(valid)
      {
	bool can_fit=(cond_numb<=cond_numb_stored);
	if(!can_fit) verbosity_lv2_master_printf(" Condition number %lg>=%lg (stored)\n",cond_numb,cond_numb_stored);
	valid&=can_fit;
      }
    
    //check that approximation it is not exagerated
    if(valid)
      {
	double exc=pow(enl_gen,4);
	double fact=cond_numb_stored/cond_numb;
	bool does_not_exceed=(fact<exc);
	if(!does_not_exceed) verbosity_lv2_master_printf(" Condition number %lg more than %lg smaller (%lg) than %lg (stored)\n",cond_numb,exc,fact,cond_numb_stored);
	valid&=does_not_exceed;
      }
    
    //check that the approximation accuracy is close to the required one
    if(valid)
      {
	double toll=2;
	bool is_not_too_precise=(appr.maxerr>=maxerr/toll);
	bool is_not_too_inaccur=(appr.maxerr<=maxerr*toll);
	bool is_accurate=(is_not_too_precise&&is_not_too_inaccur);
	if(!is_accurate) verbosity_lv2_master_printf(" Accuracy %lg more than %lg larger or smaller (%lg) than  %lg (stored)\n",
						     maxerr,toll,appr.maxerr/maxerr,appr.maxerr);
	valid&=is_accurate;
      }
    
    return valid;
  }
  
  //scale the rational expansion
  THREADABLE_FUNCTION_4ARG(set_expansions, std::vector<rat_approx_t>*,rat_appr, quad_su3**,eo_conf, theory_pars_t*,theory_pars, hmc_evol_pars_t*,evol_pars)
  {
    GET_THREAD_ID();
    
    //loop over each flav
    int nflavs=theory_pars->nflavs();
    
    //list of rat_approx to recreate
    int nto_recreate=0;
    int iappr_to_recreate[nappr_per_quark*nflavs];
    double min_to_recreate[nappr_per_quark*nflavs];
    double max_to_recreate[nappr_per_quark*nflavs];
    double maxerr_to_recreate[nappr_per_quark*nflavs];
    
    //allocate or not clover term and inverse evn clover term
    clover_term_t *Cl[2]={NULL,NULL};
    if(theory_pars->clover_to_be_computed())
      {
	for(int eo=0;eo<2;eo++) Cl[eo]=nissa_malloc("Cl",loc_volh,clover_term_t);
	chromo_operator(Cl,eo_conf);
      }
    
    //check that we have the appropriate number of quarks
    THREAD_ATOMIC_EXEC(if(IS_MASTER_THREAD) rat_appr->resize(nappr_per_quark*nflavs));
    
    const int max_iter=1000;
    for(int iflav=0;iflav<nflavs;iflav++)
      {
	quark_content_t &q=theory_pars->quarks[iflav];
	
	//find min and max eigenvalue
	double eig_min,eig_max;
	max_eigenval(&eig_max,&q,eo_conf,Cl,theory_pars->backfield[iflav],max_iter);
	switch(q.discretiz)
	  {
	  case ferm_discretiz::ROOT_STAG:
	    eig_min=pow(q.mass,2);
	    break;
	  case ferm_discretiz::ROOT_TM_CLOV:
	    eig_min=pow(q.mass,2); //possibly understimate, that's why this is separate
	    break;
	  default: crash("unknown");
	    eig_min=0;
	  }
	
	//take the pointer to the rational approximations for current flavor and mark down degeneracy
	rat_approx_t *appr=&(*rat_appr)[nappr_per_quark*iflav];
	int deg=q.deg;
	int npf=evol_pars->npseudo_fs[iflav];
	
	//generate the three approximations
	int root=ferm_discretiz::root_needed(q.discretiz);
	int extra_fact[nappr_per_quark]={2*root,-root,-root};
	double maxerr[nappr_per_quark]={sqrt(evol_pars->pf_action_residue),sqrt(evol_pars->pf_action_residue),sqrt(evol_pars->md_residue)};
	for(int i=0;i<nappr_per_quark;i++)
	  {
	    //compute condition number and exponent
	    int num=deg,den=extra_fact[i]*npf;
	    double expo=(double)num/den;
	    
	    //print the name
	    if(IS_MASTER_THREAD) snprintf(appr[i].name,20,"x^(%d/%d)",num,den);
	    THREAD_BARRIER();
	    
	    //check if the approx is valid or not and in any case fit exponents
	    verbosity_lv2_master_printf("Checking validity of approx %d/3 for flav %d/%d (%s)\n",i+1,iflav+1,nflavs,appr[i].name);
	    bool valid=check_approx_validity(appr[i],eig_min,eig_max,expo,maxerr[i]);
	    THREAD_BARRIER();
	    appr[i].num=num;
	    appr[i].den=den;
	    THREAD_BARRIER();
	    
	    //if the approximation is valid scale it
	    if(valid)
	      {
		if(num==-den) verbosity_lv2_master_printf("Exact approximation, not modifying\n");
		else
		  {
		    verbosity_lv2_master_printf("Stored rational approximation valid, adapting it quickly\n");
		    
		    if(IS_MASTER_THREAD)
		      {
			//center the arithmetic average
			double scale=sqrt(eig_max*eig_min)/sqrt(appr[i].minimum*appr[i].maximum);
			double scale_extra=pow(scale,expo);
			
			//scale the rational approximation
			appr[i].minimum*=scale;
			appr[i].maximum*=scale;
			appr[i].cons*=scale_extra;
			for(int iterm=0;iterm<appr[i].degree();iterm++)
			  {
			    appr[i].poles[iterm]*=scale;
			    appr[i].weights[iterm]*=scale*scale_extra;
			  }
		      }
		    THREAD_BARRIER();
		  }
	      }
	    else
	      {
		verbosity_lv2_master_printf("Stored rational approximation not valid, scheduling to generate a new one\n");
		iappr_to_recreate[nto_recreate]=i+nappr_per_quark*iflav;
		min_to_recreate[nto_recreate]=eig_min/enl_gen;
		max_to_recreate[nto_recreate]=eig_max*enl_gen;
		maxerr_to_recreate[nto_recreate]=maxerr[i];
		nto_recreate++;
	      }
	  }
      }
    
    //find out who recreates what
    int rank_recreating[nto_recreate];
    int nrecreated_per_rank=(nto_recreate+nranks-1)/nranks;
    verbosity_lv1_master_printf("Need to recreate %d expansions, %d ranks available\n",nto_recreate,nranks);
    for(int ito=0;ito<nto_recreate;ito++)
      {
	rank_recreating[ito]=ito/nrecreated_per_rank;
	verbosity_lv2_master_printf(" %d will be recreated by %d\n",ito,rank_recreating[ito]);
      }
    
    //create missing ones
    for(int ito=0;ito<nto_recreate;ito++)
      if(rank_recreating[ito]==rank)
	{
	  if(IS_MASTER_THREAD and VERBOSITY_LV2) printf("Rank %d recreating approx %d\n",rank,ito);
	  rat_approx_t *rat=&(*rat_appr)[iappr_to_recreate[ito]];
	  generate_approx_of_maxerr(*rat,min_to_recreate[ito],max_to_recreate[ito],maxerr_to_recreate[ito],rat->num,rat->den);
	}
    if(IS_MASTER_THREAD) MPI_Barrier(MPI_COMM_WORLD);
    THREAD_BARRIER();
    
    //now collect from other nodes
    for(int ito=0;ito<nto_recreate;ito++)
      {
	rat_approx_t *rat=&(*rat_appr)[iappr_to_recreate[ito]];
	broadcast(rat,rank_recreating[ito]);
	verbosity_lv1_master_printf("Approximation x^(%d/%d) recreated, now %d terms present\n",rat->num,rat->den,rat->degree());
      }
    
    //wait
    if(IS_MASTER_THREAD) MPI_Barrier(MPI_COMM_WORLD);
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
}
