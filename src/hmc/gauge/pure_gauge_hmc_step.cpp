#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "hmc/gauge/gluonic_action.hpp"
#include "hmc/momenta/momenta_action.hpp"
#include "hmc/momenta/momenta_generation.hpp"
#include "hmc/gauge/pure_gauge_omelyan_integrator.hpp"
#include "new_types/new_types_definitions.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //perform a full hmc step and return the difference between final and original action
  double pure_gauge_hmc_step(quad_su3 *out_conf,quad_su3 *in_conf,theory_pars_t &theory_pars,
				 pure_gauge_evol_pars_t &evol_pars,int itraj)
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
    
    //create the momenta
    generate_hmc_momenta(H);
    
    //compute initial action
    double init_action_H=momenta_action(H);
    verbosity_lv2_master_printf("Momenta action: %lg\n",init_action_H);
    double init_action_G;
    gluonic_action(&init_action_G,out_conf,&theory_pars);
    verbosity_lv2_master_printf("Gauge action: %lg\n",init_action_G);
    double init_action=init_action_G+init_action_H;
    verbosity_lv2_master_printf("Init action: %lg\n",init_action);

    //evolve forward
    omelyan_pure_gauge_evolver(H,out_conf,&theory_pars,&evol_pars);
    
    //compute final action
    double final_action_H=momenta_action(H);
    verbosity_lv2_master_printf("Momenta action: %lg\n",final_action_H);
    double final_action_G;
    gluonic_action(&final_action_G,out_conf,&theory_pars);
    verbosity_lv2_master_printf("Gauge action: %lg\n",final_action_G);
    double final_action=final_action_G+final_action_H;
    verbosity_lv2_master_printf("Final action: %lg\n",final_action);
    
    //compute the diff
    double diff_action=final_action-init_action;
    
    //free stuff
    nissa_free(H);
    
    //take time
    hmc_time+=take_time();
    verbosity_lv1_master_printf("Total time to perform pure gauge hmc step: %lg s\n",hmc_time);
    
    return diff_action;
  }
}
