#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <math.h>
#include <vector>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/rat_approx.hpp"
#include "operations/smearing/stout.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"

#include "hmc/momenta/momenta_generation.hpp"
#include "hmc/backfield.hpp"

#include "rootst_eoimpr_action.hpp"
#include "rootst_eoimpr_eigenvalues.hpp"
#include "rootst_eoimpr_omelyan_integrator.hpp"
#include "rootst_eoimpr_pseudofermions_generation.hpp"
#include "rat_expansion_database.cpp"

namespace nissa
{
  //perform a full hmc step and return the difference between final and original action
  double rootst_eoimpr_rhmc_step(quad_su3 **out_conf,quad_su3 **in_conf,theory_pars_t &theory_pars,
				 hmc_evol_pars_t &simul_pars,int itraj)
  {
    //header
    master_printf("Trajectory %d (nmd: %d, ngss: %d)\n",itraj,simul_pars.nmd_steps,simul_pars.ngauge_substeps);
    master_printf("-------------------------------\n");
    
    //take initial time
    double hmc_time=-take_time();
    
    //allocate the momenta
    quad_su3 *H[2];
    for(int par=0;par<2;par++)
      H[par]=nissa_malloc("H",loc_volh,quad_su3);
    
    //B momenta
    double *H_B=NULL;
    if(theory_pars.em_field_pars.flag==2 && itraj>theory_pars.em_field_pars.meta_skip) H_B=nissa_malloc("H_B",1,double);
    
    //copy the old conf into the new
    for(int par=0;par<2;par++)
      {
	memcpy(out_conf[par],in_conf[par],loc_volh*sizeof(quad_su3));
	set_borders_invalid(out_conf[par]);
      }
    
    //initialize rational approximation for pf/action and force calculation
    rat_approx_t rat_exp_pfgen[theory_pars.nflavs],rat_exp_actio[theory_pars.nflavs];
    for(int iflav=0;iflav<theory_pars.nflavs;iflav++)
      {
	rat_approx_create(&(rat_exp_pfgen[iflav]),db_rat_exp_nterms,"pfgen");
	rat_approx_create(&(rat_exp_actio[iflav]),db_rat_exp_nterms,"actio");
      }
    
    //allocate pseudo-fermions
    color **pf=nissa_malloc("pf*",theory_pars.nflavs,color*);
    for(int iflav=0;iflav<theory_pars.nflavs;iflav++) pf[iflav]=nissa_malloc("pf",loc_volh,color);
    
    //if needed smear the configuration for pseudo-fermions, approx generation and action computation
    //otherwise bind out_conf to sme_conf
    quad_su3 *sme_conf[2];
    for(int eo=0;eo<2;eo++) sme_conf[eo]=(theory_pars.stout_pars.nlev!=0)?
			      nissa_malloc("sme_conf",loc_volh+bord_volh+edge_volh,quad_su3):out_conf[eo];
    if(theory_pars.stout_pars.nlev!=0)
      {
	verbosity_lv2_master_printf("Stouting the links for pseudo-fermions generation and initial action computation\n");
	stout_smear(sme_conf,out_conf,&(theory_pars.stout_pars));
	
	verbosity_lv2_master_printf("Original plaquette: %16.16lg\n",global_plaquette_eo_conf(out_conf));
	verbosity_lv2_master_printf("Stouted plaquette: %16.16lg\n",global_plaquette_eo_conf(sme_conf));
	
	addrem_stagphases_to_eo_conf(sme_conf);
      }
    
    //add the phases to the output conf
    addrem_stagphases_to_eo_conf(out_conf);
    
    //generate the appropriate expansion of rational approximations
    rootst_eoimpr_scale_expansions(rat_exp_pfgen,rat_exp_actio,sme_conf,&theory_pars);
    
    //create pseudo-fermions
    for(int iflav=0;iflav<theory_pars.nflavs;iflav++)
      generate_pseudo_fermion(pf[iflav],sme_conf,theory_pars.backfield[iflav],&(rat_exp_pfgen[iflav]),
			      simul_pars.pf_action_residue);
    
    //create the momenta
    generate_hmc_momenta(H);
    if(H_B) generate_hmc_B_momenta(H_B);
    
    //compute initial action
    double init_action;
    full_rootst_eoimpr_action(&init_action,out_conf,sme_conf,H,H_B,pf,
			      &theory_pars,rat_exp_actio,simul_pars.pf_action_residue);
    
    //evolve forward
    omelyan_rootst_eoimpr_evolver(H,H_B,out_conf,pf,&theory_pars,rat_exp_actio,&simul_pars);
    
    //if needed, resmear the conf, otherwise sme_conf is already binded to out_conf
    if(theory_pars.stout_pars.nlev!=0)
      {
	verbosity_lv2_master_printf("Stouting the links for final action computation\n");
	addrem_stagphases_to_eo_conf(out_conf);
	stout_smear(sme_conf,out_conf,&(theory_pars.stout_pars));
	addrem_stagphases_to_eo_conf(out_conf);
	addrem_stagphases_to_eo_conf(sme_conf);
      }
    
    //compute final action using sme_conf (see previous note)
    double final_action;
    full_rootst_eoimpr_action(&final_action,out_conf,sme_conf,H,H_B,pf,
			      &theory_pars,rat_exp_actio,simul_pars.pf_action_residue);
    verbosity_lv2_master_printf("Final action: %lg\n",final_action);
    
    //compute the diff
    double diff_action=final_action-init_action;
    
    //remove the phases
    addrem_stagphases_to_eo_conf(out_conf);
    
    //free stuff
    for(int iflav=0;iflav<theory_pars.nflavs;iflav++) nissa_free(pf[iflav]);
    for(int par=0;par<2;par++) nissa_free(H[par]);
    if(theory_pars.stout_pars.nlev!=0) for(int eo=0;eo<2;eo++) nissa_free(sme_conf[eo]);
    nissa_free(pf);
    if(H_B) nissa_free(H_B);

    //take time
    hmc_time+=take_time();
    verbosity_lv1_master_printf("Total time to perform rhmc step: %lg s\n",hmc_time);
    
    return diff_action;
  }
}
