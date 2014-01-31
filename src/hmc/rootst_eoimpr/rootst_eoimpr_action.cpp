#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/backfield.hpp"
#include "hmc/momenta/momenta_action.hpp"
#include "hmc/gauge/tree_level_Symanzik_action.hpp"
#include "inverters/staggered/cgm_invert_stD2ee_m2.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "operations/smearing/stout.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "operations/su3_paths/topological_charge.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif


namespace nissa
{
  //compute quark action for a set of quark
  THREADABLE_FUNCTION_7ARG(rootst_eoimpr_quark_action, double*,glb_action, quad_su3**,eo_conf, int,nfl, quad_u1***,u1b, color**,pf, rat_approx_t*,appr, double,residue)
  {
    //allocate chi
    color *chi_e=nissa_malloc("chi_e",loc_volh,color);
    
    //quark action
    (*glb_action)=0;
    for(int ifl=0;ifl<nfl;ifl++)
      {
	//compute chi with background field
	add_backfield_to_conf(eo_conf,u1b[ifl]);
	summ_src_and_all_inv_stD2ee_m2_cgm(chi_e,eo_conf,appr+ifl,1000000,residue,pf[ifl]);
	rem_backfield_from_conf(eo_conf,u1b[ifl]);
	
	//compute scalar product
	double flav_action;
	double_vector_glb_scalar_prod(&flav_action,(double*)chi_e,(double*)(pf[ifl]),loc_volh*6);
	(*glb_action)+=flav_action;
      }
    
    //free
    nissa_free(chi_e);
  }
  THREADABLE_FUNCTION_END
  
  //Compute the topological action
  double topotential_action(quad_su3 **ext_conf,topotential_pars_t &pars)
  {
    quad_su3 *in_conf[2];
    if(pars.stout_pars.nlev==0)
      {
	in_conf[0]=ext_conf[0];
	in_conf[1]=ext_conf[1];
      }
    else
      {
	in_conf[0]=nissa_malloc("stout_conf",loc_volh+bord_volh+edge_volh,quad_su3);
	in_conf[1]=nissa_malloc("stout_conf",loc_volh+bord_volh+edge_volh,quad_su3);
	stout_smear(in_conf,ext_conf,&(pars.stout_pars));
      }

    //compute topocharge                                                                                              
    double charge;
    total_topological_charge_eo_conf(&charge,in_conf);
    
    //compute according to flag
    double topo_action=0;
    switch(pars.flag)
      {
      case 1: topo_action=charge*pars.theta;break;
      case 2:
	topo_action=0;
	for(std::vector<double>::iterator it=pars.past_values.begin();it!=pars.past_values.end();it++)
	  {
	    double Q=*it;
	    double diff=Q-charge,f=diff/pars.width,cont=exp(-f*f/2);
	    topo_action+=cont;
	    master_printf("Contribution: %lg\n",cont); //debug
	  }
	topo_action*=pars.coeff;
	break;
      default: crash("unknown flag %d",pars.flag);
      }
    
    return topo_action;
  }
  
  //Compute the total action of the rooted staggered e/o improved theory.
  //Passed conf must NOT contain the backfield, but contains the stagphases so remove it.
  THREADABLE_FUNCTION_8ARG(full_rootst_eoimpr_action, double*,tot_action, quad_su3**,eo_conf, quad_su3**,sme_conf, quad_su3**,H, color**,pf, theory_pars_t*,theory_pars, rat_approx_t*,appr, double,residue)
  {
    verbosity_lv1_master_printf("Computing action\n");
    
    //compute the three parts of the action
    double quark_action;
    rootst_eoimpr_quark_action(&quark_action,sme_conf,theory_pars->nflavs,theory_pars->backfield,pf,appr,residue);
    verbosity_lv1_master_printf("Quark_action: %16.16lg\n",quark_action);
    
    //gauge action
    double gluon_action;
    switch(theory_pars->gauge_action_name)
      {
      case WILSON_GAUGE_ACTION:gluon_action=theory_pars->beta*6*(1+global_plaquette_eo_conf(eo_conf))*glb_vol;break;
      case TLSYM_GAUGE_ACTION:tree_level_Symanzik_action(&gluon_action,eo_conf,theory_pars->beta,1);break;
      default:crash("Unknown action");
      }
    
    verbosity_lv1_master_printf("Gluon_action: %16.16lg\n",gluon_action);
    
    //momenta action
    double mom_action=momenta_action(H);
    verbosity_lv1_master_printf("Mom_action: %16.16lg\n",mom_action);
    
    //compute the topological action, if needed
    double topo_action=(theory_pars->topotential_pars.flag?topotential_action(eo_conf,theory_pars->topotential_pars):0);
    if(theory_pars->topotential_pars.flag) verbosity_lv1_master_printf("Topological_action: %16.16lg\n",topo_action);
    
    (*tot_action)=quark_action+gluon_action+mom_action+topo_action;
  }
  THREADABLE_FUNCTION_END
}
