#include <string.h>

#include "../../base/global_variables.h"
#include "../../base/routines.h"
#include "../../geometry/geometry_eo.h"
#include "../../base/vectors.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/rat_exp.h"

#include "../momenta/momenta_generation.h"

#include "rootst_eoimpr_action.h"
#include "rootst_eoimpr_eigenvalues.h"
#include "rootst_eoimpr_omelyan_integrator.h"
#include "rootst_eoimpr_pseudofermions_generation.h"
#include "rat_expansion_database.cpp"

//perform a full hmc step and return the difference between final and original action
double rootst_eoimpr_rhmc_step(quad_su3 **out_conf,quad_su3 **in_conf,theory_pars *physic,hmc_evol_pars *simul)
{
  double start_time=take_time();
  
  //allocate the momenta
  quad_su3 *H[2];
  for(int par=0;par<2;par++)
    H[par]=nissa_malloc("H",loc_volh,quad_su3);
  
  //copy the old conf into the new
  for(int par=0;par<2;par++)
    {
      memcpy(out_conf[par],in_conf[par],loc_volh*sizeof(quad_su3));
      set_borders_invalid(out_conf[par]);
    }
  
  //initialize rational approximation for pf/action and force calculation
  rat_approx rat_exp_pfgen[physic->nflavs],rat_exp_actio[physic->nflavs];
  for(int iflav=0;iflav<physic->nflavs;iflav++)
    {
      rat_approx_create(&(rat_exp_pfgen[iflav]),db_rat_exp_nterms,"pfgen");
      rat_approx_create(&(rat_exp_actio[iflav]),db_rat_exp_nterms,"actio");
    }

  //allocate pseudo-fermions
  color **pf=nissa_malloc("pf*",physic->nflavs,color*);
  for(int iflav=0;iflav<physic->nflavs;iflav++) pf[iflav]=nissa_malloc("pf",loc_volh,color);
  
  //add the phases
  addrem_stagphases_to_eo_conf(out_conf);
  
  //generate the appropriate expansion of rational approximations
  rootst_eoimpr_scale_expansions(rat_exp_pfgen,rat_exp_actio,out_conf,physic);
  
  //create the momenta
  generate_hmc_momenta(H);
  
  //create pseudo-fermions
  for(int iflav=0;iflav<physic->nflavs;iflav++)
    generate_pseudo_fermion(pf[iflav],out_conf,physic->backfield[iflav],&(rat_exp_pfgen[iflav]),simul->pf_action_residue);
  
  //compute initial action
  double init_action=full_rootst_eoimpr_action(out_conf,H,pf,physic,rat_exp_actio,simul->pf_action_residue);
  verbosity_lv2_master_printf("Init action: %lg\n",init_action);
  
  //evolve forward
  omelyan_rootst_eoimpr_evolver(H,out_conf,pf,physic,rat_exp_actio,simul);
  
  //compute final action
  double final_action=full_rootst_eoimpr_action(out_conf,H,pf,physic,rat_exp_actio,simul->pf_action_residue);
  verbosity_lv2_master_printf("Final action: %lg\n",final_action);
  
  //compute the diff
  double diff_action=final_action-init_action;

  //remove the phases
  addrem_stagphases_to_eo_conf(out_conf);
  
  //free stuff
  for(int iflav=0;iflav<physic->nflavs;iflav++) nissa_free(pf[iflav]);
  for(int par=0;par<2;par++) nissa_free(H[par]);
  nissa_free(pf);

  verbosity_lv1_master_printf("Total time to perform rhmc step: %lg s\n",take_time()-start_time);
  
  return diff_action;
}
