#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "hmc/gauge/gluonic_action.hpp"
#include "hmc/gauge/MFACC_fields.hpp"
#include "hmc/gauge/MFACC_fields.hpp"
#include "hmc/momenta/momenta_action.hpp"
#include "hmc/momenta/momenta_generation.hpp"
#include "hmc/gauge/pure_gauge_Omelyan_integrator.hpp"
#include "new_types/rat_approx.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  double pure_gauge_action(quad_su3 *conf,theory_pars_t &theory_pars,pure_gauge_evol_pars_t &evol_pars,quad_su3 *H,su3 **phi,su3 **pi)
  {
    //compute action for momenta
    double action_H=0;
    if(evol_pars.use_Facc) action_H=momenta_action_with_FACC(conf,evol_pars.kappa,100000,evol_pars.residue,H);
    else action_H=momenta_action(H);
    verbosity_lv1_master_printf("Momenta action: %lg\n",action_H);
    
    //compute action for FACC
    double action_phi=0,action_pi=0;
    if(evol_pars.use_Facc)
      {
	action_phi=MFACC_fields_action(phi);
	verbosity_lv1_master_printf("Fourier acceleration fields action: %lg\n",action_phi);
	action_pi=MFACC_momenta_action(pi,conf,evol_pars.kappa);
	verbosity_lv1_master_printf("Fourier acceleration momenta action: %lg\n",action_pi);
      }
    
    //compute action for G
    double action_G=0;
    gluonic_action(&action_G,conf,theory_pars.gauge_action_name,theory_pars.beta);
    verbosity_lv1_master_printf("Gauge action: %lg\n",action_G);
    
    return action_G+action_H+action_phi+action_pi;
  }
  
  //perform a full hmc step and return the difference between final and original action
  double pure_gauge_hmc_step(quad_su3 *out_conf,quad_su3 *in_conf,theory_pars_t &theory_pars,pure_gauge_evol_pars_t &evol_pars,rat_approx_t *rat_exp_H,int itraj)
  {
    if(in_conf==out_conf) crash("in==out");
    
    static int init=1;
    
    //header
    double hmc_time=-take_time();
    master_printf("Trajectory %d (nmd: %d)\n",itraj,evol_pars.nmd_steps);
    
    //allocate the momenta and copy old conf into new
    quad_su3 *H=nissa_malloc("H",loc_vol,quad_su3);
    vector_copy(out_conf,in_conf);
    
    //if we accelerate draw also momenta and position
    static su3 *phi[n_FACC_fields],*pi[n_FACC_fields];
    if(evol_pars.use_Facc)
      {
	
	if(init)
	  {
	    init=0;
	    //if we need to accelerate allocate auxiliary fields
	    for(int id=0;id<n_FACC_fields;id++)
	      {
		phi[id]=nissa_malloc("phi",loc_vol+bord_vol,su3);
		pi[id]=nissa_malloc("pi",loc_vol+bord_vol,su3);
	      }
	    
	    //generate FACC fields and momenta
	    for(int id=0;id<n_FACC_fields;id++) generate_MFACC_fields(phi[id]);
	    generate_MFACC_momenta(pi,out_conf,evol_pars.kappa,evol_pars.residue);
	  }
	
	//create the momenta
	generate_hmc_momenta_with_FACC(H,in_conf,rat_exp_H,evol_pars.kappa,evol_pars.residue);
      }
    else generate_hmc_momenta(H);
    
    //compute the action
    double init_action=pure_gauge_action(out_conf,theory_pars,evol_pars,H,phi,pi);
    verbosity_lv2_master_printf("Init action: %lg\n",init_action);
    
    //evolve forward
    if(evol_pars.use_Facc) Omelyan_pure_gauge_FACC_evolver(H,out_conf,pi,phi,&theory_pars,&evol_pars);
    else                   Omelyan_pure_gauge_evolver(H,out_conf,&theory_pars,&evol_pars);
    
    //compute the action
    double final_action=pure_gauge_action(out_conf,theory_pars,evol_pars,H,phi,pi);
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
      for(int id=0;id<n_FACC_fields;id++)
	{
	  //nissa_free(phi[id]);
	  //nissa_free(pi[id]);
	}
    
    return diff_action;
  }
}
