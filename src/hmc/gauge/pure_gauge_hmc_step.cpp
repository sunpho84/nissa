#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "hmc/gauge/gluonic_action.hpp"
#include "hmc/gauge/MFACC_fields.hpp"
#include "hmc/gauge/MFACC_fields.hpp"
#include "hmc/momenta/momenta_action.hpp"
#include "hmc/momenta/momenta_generation.hpp"
#include "hmc/gauge/pure_gauge_omelyan_integrator.hpp"
#include "new_types/new_types_definitions.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //perform a full hmc step and return the difference between final and original action
  double pure_gauge_hmc_step(quad_su3 *out_conf,quad_su3 *in_conf,theory_pars_t &theory_pars,pure_gauge_evol_pars_t &evol_pars,rat_approx_t *rat_exp_H,int itraj)
  {
    if(in_conf==out_conf) crash("in==out");
    
    //header
    master_printf("Trajectory %d (nmd: %d)\n",itraj,evol_pars.nmd_steps);
    
    //take initial time
    double hmc_time=-take_time();
    
    //allocate the momenta
    quad_su3 *H=nissa_malloc("H",loc_vol,quad_su3);
    
    //copy the old conf into the new
    vector_copy(out_conf,in_conf);
    
    //if we accelerate draw also momenta and position
    su3 *phi[2],*pi[2];
    if(evol_pars.use_Facc)
      {
	//if we need to accelerate allocate auxiliary fields
	if(evol_pars.use_Facc)
	  for(int id=0;id<2;id++)
	    {
	      phi[id]=nissa_malloc("phi",loc_vol+bord_vol,su3);
	      pi[id]=nissa_malloc("pi",loc_vol+bord_vol,su3);
	    }
	
	//create the momenta
	generate_hmc_momenta_with_FACC(H,in_conf,rat_exp_H,evol_pars.kappa,evol_pars.residue);
	
	for(int id=0;id<2;id++) generate_MFACC_fields(phi[id]);
	generate_MFACC_momenta(pi,out_conf,evol_pars.kappa,evol_pars.residue);
      }
    else generate_hmc_momenta(H);

    
    //compute initial action
    double init_action_H=momenta_action(H);
    verbosity_lv2_master_printf("Momenta action: %lg\n",init_action_H);
    double init_action_phi=0,init_action_pi=0;
    if(evol_pars.use_Facc)
      {
	init_action_phi=MFACC_fields_action(phi);
	verbosity_lv2_master_printf("Fourier acceleration fields action: %lg\n",init_action_phi);
	init_action_pi=MFACC_momenta_action(pi,out_conf,evol_pars.kappa);
	verbosity_lv2_master_printf("Fourier acceleration momenta action: %lg\n",init_action_pi);
      }
    double init_action_G;
    gluonic_action(&init_action_G,out_conf,&theory_pars,false);
    verbosity_lv2_master_printf("Gauge action: %lg\n",init_action_G);
    double init_action=init_action_G+init_action_H+init_action_phi+init_action_pi;
    verbosity_lv2_master_printf("Init action: %lg\n",init_action);
    
    //evolve forward
    if(evol_pars.use_Facc) omelyan_pure_gauge_Facc_evolver(H,out_conf,pi,phi,&theory_pars,&evol_pars);
    else                   omelyan_pure_gauge_evolver(H,out_conf,&theory_pars,&evol_pars);
    
    //compute final action
    double final_action_H=momenta_action(H);
    verbosity_lv2_master_printf("Momenta action: %lg\n",final_action_H);
    double final_action_phi=0,final_action_pi=0;
    if(evol_pars.use_Facc)
      {
	final_action_phi=MFACC_fields_action(phi);
	verbosity_lv2_master_printf("Fourier acceleration fields action: %lg\n",final_action_phi);
	final_action_pi=MFACC_momenta_action(pi,out_conf,evol_pars.kappa);
	verbosity_lv2_master_printf("Fourier acceleration momenta action: %lg\n",final_action_pi);
      }
    double final_action_G;
    gluonic_action(&final_action_G,out_conf,&theory_pars,false);
    verbosity_lv2_master_printf("Gauge action: %lg\n",final_action_G);
    double final_action=final_action_G+final_action_H+final_action_phi+final_action_pi;
    verbosity_lv2_master_printf("Final action: %lg\n",final_action);
    
    //compute the diff
    double diff_action=final_action-init_action;
    
    //free stuff
    nissa_free(H);
    
    //take time
    hmc_time+=take_time();
    verbosity_lv1_master_printf("Total time to perform pure gauge hmc step: %lg s\n",hmc_time);
    
    //if accelerated, free everything
    if(evol_pars.use_Facc)
      for(int id=0;id<2;id++)
	{
	  nissa_free(phi[id]);
	  nissa_free(pi[id]);
	}
    
    return diff_action;
  }
}
