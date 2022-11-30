#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "hmc/gauge/gluonic_action.hpp"
#include "hmc/momenta/momenta_action.hpp"
#include "hmc/momenta/momenta_generation.hpp"
#include "hmc/gauge/pure_gauge_Omelyan_integrator.hpp"
#include "new_types/rat_approx.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  
  double pure_gauge_action(quad_su3 *conf,theory_pars_t &theory_pars,pure_gauge_evol_pars_t &evol_pars,quad_su3 *H)
  {
    //compute action for momenta
    double action_H=momenta_action(H);
    verbosity_lv1_master_printf("Momenta action: %lg\n",action_H);
    
    //compute action for G
    double action_G=0;
    gluonic_action(&action_G,conf,theory_pars.gauge_action_name,theory_pars.beta);
    verbosity_lv1_master_printf("Gauge action: %lg\n",action_G);
    
    return action_G+action_H;
  }
  
  //perform a full hmc step and return the difference between final and original action
  double pure_gauge_hmc_step(quad_su3 *out_conf,quad_su3 *in_conf,theory_pars_t &theory_pars,pure_gauge_evol_pars_t &evol_pars,rat_approx_t *rat_exp_H,int itraj)
  {
    if(in_conf==out_conf) crash("in==out");
    
    //header
    double hmc_time=-take_time();
    master_printf("Trajectory %d (nmd: %d)\n",itraj,evol_pars.nmd_steps);
    
    //allocate the momenta and copy old conf into new
    quad_su3 *H=nissa_malloc("H",locVol,quad_su3);
    vector_copy(out_conf,in_conf);
    
    generate_hmc_momenta(H);
    
    //compute the action
    double init_action=pure_gauge_action(out_conf,theory_pars,evol_pars,H);
    verbosity_lv2_master_printf("Init action: %lg\n",init_action);
    
    Omelyan_pure_gauge_evolver(H,out_conf,&theory_pars,&evol_pars);
    
    //compute the action
    double final_action=pure_gauge_action(out_conf,theory_pars,evol_pars,H);
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
