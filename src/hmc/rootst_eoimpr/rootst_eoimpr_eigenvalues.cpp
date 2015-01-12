#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "communicate/communicate.hpp"
#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "inverters/staggered/cg_invert_stD2ee_m2.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/rat_approx.hpp"
#include "new_types/su3.hpp"
#include "operations/remez/remez_algorithm.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "hmc/backfield.hpp"

#define ENL_GEN 1.18920711500272106671

namespace nissa
{
  //Return the maximal eigenvalue of the staggered Dirac operator for the passed quark
  //assumes that the passed conf already has stag phases inside it
  THREADABLE_FUNCTION_4ARG(max_eigenval, double*,eig_max, quark_content_t*,quark_content, quad_su3**,eo_conf, int,niters)
  {
    niters=50;
    GET_THREAD_ID();
    communicate_ev_and_od_quad_su3_borders(eo_conf);
    (*eig_max)=0;
    
    color *vec_in=nissa_malloc("vec_in",loc_volh+bord_volh,color);
    color *vec_out=nissa_malloc("vec_out",loc_volh,color);
    color *tmp=nissa_malloc("tmp",loc_volh+bord_volh,color);
    
    //generate the random field
    if(IS_MASTER_THREAD)
      NISSA_LOC_VOLH_LOOP(ivol)
	color_put_to_gauss(vec_in[ivol],&(loc_rnd_gen[loclx_of_loceo[EVN][ivol]]),3);
    
    set_borders_invalid(vec_in);
    
    //perform initial normalization
    double init_norm;
    double_vector_normalize(&init_norm,(double*)vec_in,(double*)vec_in,glb_volh*3,6*loc_volh);
    verbosity_lv3_master_printf("Init norm: %lg\n",init_norm);
    
    //apply the vector niter times normalizing at each iter
    for(int iter=0;iter<niters;iter++)
      {
	//inv_stD2ee_m2_cg(vec_out,NULL,eo_conf,sqr(quark_content.mass),10000,5,1.e-13,vec_in);
	apply_stD2ee_m2(vec_out,eo_conf,tmp,quark_content->mass*quark_content->mass,vec_in);
	
	//compute the norm
	double_vector_normalize(eig_max,(double*)vec_in,(double*)vec_out,glb_volh*3,6*loc_volh);
	
	verbosity_lv2_master_printf("max_eigen search mass %lg, iter=%d, eig=%16.16lg\n",quark_content->mass,iter,*eig_max);
      }
    
    //assume a 10% excess
    (*eig_max)*=1.1;
    
    nissa_free(tmp);
    nissa_free(vec_out);
    nissa_free(vec_in);
    
    verbosity_lv2_master_printf("max_eigen mass %lg: %16.16lg\n",quark_content->mass,*eig_max);
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
	bool allocated=(appr.degree!=0);
	if(!allocated) verbosity_lv3_master_printf(" Not allocated\n");
	valid&=allocated;
      }
    
    //check that exponent is the same
    if(valid)
      {
	double expo_stored=(double)appr.num/appr.den;
	bool expo_is_the_same=(fabs(expo-expo_stored)<1e-14);
	if(!expo_is_the_same) verbosity_lv3_master_printf(" Approx stored is for x^%lg, required is for x^%lg\n",expo_stored,expo);
	valid&=expo_is_the_same;
      }
    
    //check that it can be fitted
    if(valid)
      {
	bool can_fit=(cond_numb<=cond_numb_stored);
	if(!can_fit) verbosity_lv3_master_printf(" Condition number %lg>=%lg (stored)\n",cond_numb,cond_numb_stored);
	valid&=can_fit;
      }
    
    //check that approximation it is not exagerated
    if(valid)
      {
	double exc=pow(ENL_GEN,4);
	double fact=cond_numb_stored/cond_numb;
	bool does_not_exceed=(fact<exc);
	if(!does_not_exceed) verbosity_lv3_master_printf(" Condition number %lg more than %lg smaller (%lg) than %lg (stored)\n",cond_numb,exc,fact,cond_numb_stored);
	valid&=does_not_exceed;
      }
    
    //check that the approximation accuracy is close to the required one
    if(valid)
      {
	double toll=2;
	bool is_not_too_precise=(appr.maxerr>=maxerr/toll);
	bool is_not_too_inaccur=(appr.maxerr<=maxerr*toll);
	bool is_accurate=(is_not_too_precise&&is_not_too_inaccur);
	if(!is_accurate) verbosity_lv3_master_printf(" Accuracy %lg more than %lg larger or smaller (%lg) than  %lg (stored)\n",
						     maxerr,toll,appr.maxerr/maxerr,appr.maxerr);
	valid&=is_accurate;
      }
    
    return valid;
  }
  
  //scale the rational expansion
  //assumes that the conf has already stag phases inside
  THREADABLE_FUNCTION_3ARG(rootst_eoimpr_set_expansions, hmc_evol_pars_t*,evol_pars, quad_su3**,eo_conf, theory_pars_t*,theory_pars)
  {
    GET_THREAD_ID();
    
    //loop over each flav
    int nflavs=theory_pars->nflavs;
    for(int iflav=0;iflav<nflavs;iflav++)
      {
	//find min and max eigenvalue
	double eig_min,eig_max;
	add_backfield_to_conf(eo_conf,theory_pars->backfield[iflav]);
	max_eigenval(&eig_max,&(theory_pars->quark_content[iflav]),eo_conf,50);
	eig_min=pow(theory_pars->quark_content[iflav].mass,2);
	rem_backfield_from_conf(eo_conf,theory_pars->backfield[iflav]);
	
	//take the pointer to the rational approximations for current flavor and mark down degeneracy
	rat_approx_t *appr=evol_pars->rat_appr+3*iflav;
	int deg=theory_pars->quark_content[iflav].deg;
	//the number of pseudofermion is set to 1 for the time being
	int npf=evol_pars->npseudo_fs[iflav]=1;

	//generate the three approximations
	int extra_fact[3]={8,-4,-4};
	double maxerr[3]={sqrt(evol_pars->pf_action_residue),sqrt(evol_pars->pf_action_residue),sqrt(evol_pars->md_residue)};
	for(int i=0;i<3;i++)
	  {
	    //compute condition number and exponent
	    int num=deg,den=extra_fact[i]*npf;
	    double expo=(double)num/den;
	    
	    //print the name
	    if(IS_MASTER_THREAD) snprintf(appr[i].name,20,"x^(%d/%d)",num,den);
	    THREAD_BARRIER();
	    
	    //check if the approx is valid or not and in any case fit exponents	    
	    verbosity_lv2_master_printf("Checking validity of approx %s (%d/3) for flav %d/%d\n",appr[i].name,i,iflav+1,nflavs);
	    bool valid=check_approx_validity(appr[i],eig_min,eig_max,expo,maxerr[i]);
	    THREAD_BARRIER();
	    appr[i].num=num;
	    appr[i].den=den;
	    THREAD_BARRIER();
	    
	    //if the approximation is valid scale it	    
	    if(valid)
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
		    for(int iterm=0;iterm<appr[i].degree;iterm++)
		      {
			appr[i].poles[iterm]*=scale;
			appr[i].weights[iterm]*=scale*scale_extra;
		      }
		  }
		THREAD_BARRIER();
	      }
	    else
	      {
		verbosity_lv2_master_printf("Stored rational approximation not valid, generating a new one\n");
		generate_approx_of_maxerr(appr[i],eig_min/ENL_GEN,eig_max*ENL_GEN,maxerr[i],num,den);
	      }
	    
	    if(VERBOSITY_LV3) master_printf_rat_approx(appr+i);
	  }
      }
  }
  THREADABLE_FUNCTION_END
  
}
