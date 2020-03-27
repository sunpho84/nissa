#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "hmc/backfield.hpp"
#include "hmc/gauge/gluonic_action.hpp"
#include "hmc/gauge/topological_action.hpp"
#include "hmc/hmc.hpp"
#include "hmc/momenta/momenta_action.hpp"
#include "hmc/multipseudo/multipseudo_rhmc_step.hpp"
#include "geometry/geometry_eo.hpp"
#include "inverters/staggered/cgm_invert_stD2ee_m2.hpp"
#include "inverters/twisted_clover/cgm_invert_tmclovDkern_eoprec_square_portable.hpp"
#include "linalgs/linalgs.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "threads/threads.hpp"

#include <algorithm>

namespace nissa
{
  //compute quark action for a set of quark
  THREADABLE_FUNCTION_7ARG(compute_quark_action, double*,glb_action, quad_su3**,eo_conf, std::vector<quad_u1**>,u1b, std::vector<std::vector<pseudofermion_t> >*,pf, std::vector<quark_content_t>,quark_content, hmc_evol_pars_t*,simul_pars, std::vector<rat_approx_t>*,rat_appr)
  {
    int nfl=quark_content.size();
    double res=simul_pars->pf_action_residue;
    
    //allocate or not clover term and inverse evn clover term
    clover_term_t *Cl[2]={NULL,NULL};
    inv_clover_term_t *invCl_evn=NULL;
    bool clover_to_be_computed=false;
    for(int iflav=0;iflav<nfl;iflav++) clover_to_be_computed|=ferm_discretiz::include_clover(quark_content[iflav].discretiz);
    if(clover_to_be_computed)
      {
	for(int eo=0;eo<2;eo++) Cl[eo]=nissa_malloc("Cl",loc_volh,clover_term_t);
	invCl_evn=nissa_malloc("invCl_evn",loc_volh,inv_clover_term_t);
	chromo_operator(Cl,eo_conf);
      }
    
    //quark action
    (*glb_action)=0;
    for(int ifl=0;ifl<nfl;ifl++)
      {
	rat_approx_t *r=&((*rat_appr)[nappr_per_quark*ifl+RAT_APPR_QUARK_ACTION]);
	quark_content_t &q=quark_content[ifl];
	pseudofermion_t chi_e(q.discretiz,"chi_e");
	
	//if clover is included, compute it
	if(clover_to_be_computed)
	  {
	    chromo_operator_include_cSW(Cl,q.cSW);
	    invert_twisted_clover_term(invCl_evn,q.mass,q.kappa,Cl[EVN]);
	  }
	
	for(int ipf=0;ipf<simul_pars->npseudo_fs[ifl];ipf++)
	  {
	    pseudofermion_t &p=(*pf)[ifl][ipf];
	    verbosity_lv1_master_printf("Computing action for flavour %d/%d, pseudofermion %d/%d\n",ifl+1,nfl,ipf+1,simul_pars->npseudo_fs[ifl]);
	    
	    //compute chi with background field
	    switch(q.discretiz)
	      {
	      case ferm_discretiz::ROOT_STAG:
		add_backfield_with_stagphases_to_conf(eo_conf,u1b[ifl]);
		summ_src_and_all_inv_stD2ee_m2_cgm(chi_e.stag,eo_conf,r,1000000,res,p.stag);
		rem_backfield_with_stagphases_from_conf(eo_conf,u1b[ifl]);
		break;
	      case ferm_discretiz::ROOT_TM_CLOV:
		add_backfield_without_stagphases_to_conf(eo_conf,u1b[ifl]);
		summ_src_and_all_inv_tmclovDkern_eoprec_square_portable(chi_e.Wils,eo_conf,q.kappa,Cl[ODD],invCl_evn,q.mass,r,1000000,res,p.Wils);
		rem_backfield_without_stagphases_from_conf(eo_conf,u1b[ifl]);
		break;
	      default: crash("still not implemented");
	      }
	    
	    //compute scalar product
	    double flav_action=chi_e.scal_prod_with(p);
	    (*glb_action)+=flav_action;
	  }
	
	//remove cSW from chromo operator
	if(clover_to_be_computed) chromo_operator_remove_cSW(Cl,q.cSW);
      }
    
    //free clover term if ever allocated
    if(clover_to_be_computed)
      {
	nissa_free(invCl_evn);
	for(int eo=0;eo<2;eo++) nissa_free(Cl[eo]);
      }
  }
  THREADABLE_FUNCTION_END
  
  //Compute the total action of the rooted staggered e/o improved theory
  THREADABLE_FUNCTION_9ARG(full_theory_action, double*,tot_action, quad_su3**,eo_conf, quad_su3**,sme_conf, quad_su3**,H, std::vector<std::vector<pseudofermion_t> > *,pf, theory_pars_t*,theory_pars, hmc_evol_pars_t*,simul_pars, std::vector<rat_approx_t>*,rat_appr, double,external_quark_action)
  {
    verbosity_lv1_master_printf("Computing action\n");
    
    //compute the three parts of the action
    double quark_action;
    if(external_quark_action>=0)
      {
	verbosity_lv1_master_printf("No need to compute pseudofermion action\n");
	quark_action=external_quark_action;
      }
    else compute_quark_action(&quark_action,sme_conf,theory_pars->backfield,pf,theory_pars->quarks,simul_pars,rat_appr);
    verbosity_lv1_master_printf("Quark_action: %16.16lg\n",quark_action);
    
    //gauge action
    double gluon_action;
    gluonic_action(&gluon_action,eo_conf,theory_pars->gauge_action_name,theory_pars->beta);
    verbosity_lv1_master_printf("Gluon_action: %16.16lg\n",gluon_action);
    
    //momenta action
    double mom_action=momenta_action(H);
    verbosity_lv1_master_printf("Mom_action: %16.16lg\n",mom_action);
    
    //topological action
    double topo_action=(theory_pars->topotential_pars.flag?topotential_action(eo_conf,theory_pars->topotential_pars):0);
    if(theory_pars->topotential_pars.flag) verbosity_lv1_master_printf("Topological_action: %16.16lg\n",topo_action);
    
    //total action
    (*tot_action)=quark_action+gluon_action+mom_action+topo_action;
  }
  THREADABLE_FUNCTION_END
}
