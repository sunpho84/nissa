#include "../../base/debug.h"
#include "../../base/global_variables.h"
#include "../../base/routines.h"
#include "../../base/vectors.h"
#include "../../geometry/geometry_eo.h"
#include "../../inverters/staggered/cgm_invert_stD2ee_m2.h"
#include "../../new_types/new_types_definitions.h"
#include "../../operations/su3_paths/plaquette.h"
#include "../../hmc/gauge/tree_level_Symanzik_action.h"

#include "../backfield.h"
#include "../momenta/momenta_action.h"

//compute quark action for a set of quark
double rootst_eoimpr_quark_action(quad_su3 **eo_conf,int nfl,quad_u1 ***u1b,color **pf,rat_approx *appr,double residue)
{  
  //allocate chi
  color *chi_e=nissa_malloc("chi_e",loc_volh,color);
  
  //quark action
  double loc_action=0;
  for(int ifl=0;ifl<nfl;ifl++)
    {
      //compute chi with background field
      add_backfield_to_conf(eo_conf,u1b[ifl]);
      summ_src_and_all_inv_stD2ee_m2_cgm(chi_e,eo_conf,appr+ifl,1000000,residue,pf[ifl]);
      rem_backfield_from_conf(eo_conf,u1b[ifl]);
      
      //compute scalar product
      nissa_loc_volh_loop(ivol)
	for(int ic=0;ic<3;ic++)
	  loc_action+=chi_e[ivol][ic][0]*pf[ifl][ivol][ic][0]+chi_e[ivol][ic][1]*pf[ifl][ivol][ic][1];
    }
  
  //global reducton
  double glb_action=glb_reduce_double(loc_action);
  
  //free
  nissa_free(chi_e);
  
  return glb_action;
}

//Compute the total action of the rooted staggered e/o improved theory.
//Passed conf must NOT contain the backfield, but contains the stagphases so remove it.
double full_rootst_eoimpr_action(quad_su3 **eo_conf,quad_su3 **sme_conf,quad_su3 **H,color **pf,theory_pars *physics,rat_approx *appr,double residue)
{
  verbosity_lv1_master_printf("Computing action\n");

  //compute the three parts of the action
  double quark_action=rootst_eoimpr_quark_action(sme_conf,physics->nflavs,physics->backfield,pf,appr,residue);
  verbosity_lv2_master_printf("Quark_action: %16.16lg\n",quark_action);
  
  //gauge action
  double gluon_action;
  switch(physics->gac_type)
    {
    case Wilson_action: 
      gluon_action=physics->beta*6*(1+global_plaquette_eo_conf(eo_conf))*glb_vol;
      break;
    case tlSym_action:
      addrem_stagphases_to_eo_conf(eo_conf);
      gluon_action=tree_level_Symanzik_action(eo_conf,physics->beta);
      addrem_stagphases_to_eo_conf(eo_conf);
      break;
    default:
      crash("Unknown action");
    }
  
  verbosity_lv2_master_printf("Gluon_action: %16.16lg\n",gluon_action);
  
  //momenta action
  double mom_action=momenta_action(H);
  verbosity_lv2_master_printf("Mom_action: %16.16lg\n",mom_action);
  
  return quark_action+gluon_action+mom_action;
}
