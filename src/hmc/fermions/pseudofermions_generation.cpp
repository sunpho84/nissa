#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/multipseudo/multipseudo_rhmc_step.hpp"
#include "inverters/staggered/cgm_invert_stD2ee_m2.hpp"
#include "inverters/twisted_clover/cgm_invert_tmclovDkern_eoprec_square_portable.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "threads/threads.hpp"
#include "hmc/backfield.hpp"
#include "hmc/hmc.hpp"

namespace nissa
{
  //generate pseudo-fermion using color vector generator
  THREADABLE_FUNCTION_9ARG(generate_pseudo_fermion, double*,action, pseudofermion_t*,pf, quad_su3**,conf, clover_term_t*,Cl_odd, inv_clover_term_t*,invCl_evn, quad_u1**,u1b, rat_approx_t*,rat, double,residue, quark_content_t,q)
  {
    //generate the random field
    pseudofermion_t pf_hb_vec(q.discretiz);
    pf_hb_vec.fill();
    
    //compute action
    (*action)=pf_hb_vec.norm2();
    
    //invert to perform hv
    add_backfield_with_stagphases_to_conf(conf,u1b);
    switch(q.discretiz)
      {
      case ferm_discretiz::ROOT_STAG:
	summ_src_and_all_inv_stD2ee_m2_cgm(pf->stag,conf,rat,10000000,residue,pf_hb_vec.stag);
	break;
      case ferm_discretiz::ROOT_TM_CLOV:
	summ_src_and_all_inv_tmclovDkern_eoprec_square_portable(pf->Wils,conf,q.kappa,Cl_odd,invCl_evn,q.mass,rat,10000000,residue,pf_hb_vec.Wils);
	break;
      default:crash("not supported");break;
      }
    rem_backfield_with_stagphases_from_conf(conf,u1b);
  }
  THREADABLE_FUNCTION_END
  
  //gemerate all pseudofermions
  double generate_pseudofermions(std::vector<std::vector<pseudofermion_t> > &pf,quad_su3 **conf,theory_pars_t &theory_pars,hmc_evol_pars_t &simul_pars,std::vector<rat_approx_t> &rat_appr)
  {
    //allocate or not clover term and inverse evn clover term
    clover_term_t *Cl[2]={NULL,NULL};
    inv_clover_term_t *invCl_evn=NULL;
    bool clover_to_be_computed=false;
    for(int iflav=0;iflav<theory_pars.nflavs();iflav++) clover_to_be_computed|=ferm_discretiz::include_clover(theory_pars.quarks[iflav].discretiz);
    if(clover_to_be_computed)
      {
	for(int eo=0;eo<2;eo++) Cl[eo]=nissa_malloc("Cl",loc_volh,clover_term_t);
	invCl_evn=nissa_malloc("invCl_evn",loc_volh,inv_clover_term_t);
	chromo_operator(Cl,conf);
      }
    
    //create pseudo-fermions and store action
    double pf_action=0;
    for(int iflav=0;iflav<theory_pars.nflavs();iflav++)
      {
	//get quark
	quark_content_t &q=theory_pars.quarks[iflav];
	
	//if clover is included, compute it
	if(ferm_discretiz::include_clover(q.discretiz))
	  {
	    chromo_operator_include_cSW(Cl,q.cSW);
	    invert_twisted_clover_term(invCl_evn,q.mass,q.kappa,Cl[EVN]);
	  }
	
	//loop over pseudofermions
	for(int ipf=0;ipf<simul_pars.npseudo_fs[iflav];ipf++)
	  {
	    verbosity_lv1_master_printf("Generating pseudofermion %d/%d for flavour %d/%d\n",ipf+1,simul_pars.npseudo_fs[iflav],iflav+1,theory_pars.nflavs());
	    double pf_action_flav;
	    if(invCl_evn) master_printf("Cl_odd: %lg invCl_evn: %lg\n",Cl[ODD][0][0][0][0][0],invCl_evn[0][0][0][0][0][0][0]);
	    generate_pseudo_fermion(&pf_action_flav,&(pf[iflav][ipf]),conf,Cl[ODD],invCl_evn,theory_pars.backfield[iflav],&rat_appr[nappr_per_quark*iflav+RAT_APPR_PF_GEN],simul_pars.pf_action_residue,q);
	    pf_action+=pf_action_flav;
	  }
	
	//remove cSW from chromo operator
	if(ferm_discretiz::include_clover(q.discretiz)) chromo_operator_remove_cSW(Cl,q.cSW);
      }
    
    //free clover term if ever allocated
    if(clover_to_be_computed)
      {
	nissa_free(invCl_evn);
	for(int eo=0;eo<2;eo++) nissa_free(Cl[eo]);
      }
    
    return pf_action;
  }
}
