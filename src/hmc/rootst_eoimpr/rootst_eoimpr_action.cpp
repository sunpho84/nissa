#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/debug.h"
#include "../../base/global_variables.h"
#include "../../base/thread_macros.h"
#include "../../base/vectors.h"
#include "../../geometry/geometry_eo.h"
#include "../../inverters/staggered/cgm_invert_stD2ee_m2.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/new_types_definitions.h"
#include "../../operations/su3_paths/plaquette.h"
#include "../../hmc/gauge/tree_level_Symanzik_action.h"
#include "../../routines/ios.h"
#include "../../routines/mpi_routines.h"
#ifdef USE_THREADS
 #include "../../routines/thread.h"
#endif

#include "../backfield.h"
#include "../momenta/momenta_action.h"

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
}}

//Compute the total action of the rooted staggered e/o improved theory.
//Passed conf must NOT contain the backfield, but contains the stagphases so remove it.
THREADABLE_FUNCTION_8ARG(full_rootst_eoimpr_action, double*,tot_action, quad_su3**,eo_conf, quad_su3**,sme_conf, quad_su3**,H, color**,pf, theory_pars_t*,theory_pars, rat_approx_t*,appr, double,residue)
{
  verbosity_lv1_master_printf("Computing action\n");

  //compute the three parts of the action
  double quark_action;
  rootst_eoimpr_quark_action(&quark_action,sme_conf,theory_pars->nflavs,theory_pars->backfield,pf,appr,residue);
  verbosity_lv2_master_printf("Quark_action: %16.16lg\n",quark_action);
  
  //gauge action
  double gluon_action;
  switch(theory_pars->gauge_action_name)
    {
    case Wilson_action:
      gluon_action=theory_pars->beta*6*(1+global_plaquette_eo_conf(eo_conf))*glb_vol;
      break;
    case tlSym_action:
      tree_level_Symanzik_action(&gluon_action,eo_conf,theory_pars->beta,1);
      break;
    default:
      crash("Unknown action");
    }
  
  verbosity_lv2_master_printf("Gluon_action: %16.16lg\n",gluon_action);
  
  //momenta action
  double mom_action=momenta_action(H);
  verbosity_lv2_master_printf("Mom_action: %16.16lg\n",mom_action);
  
  (*tot_action)=quark_action+gluon_action+mom_action;
}}
