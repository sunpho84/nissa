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
    CRASH("reimplement");
    
    // //compute action for momenta
    // const double action_H=momenta_action(H);
    // VERBOSITY_LV1_MASTER_PRINTF("Momenta action: %lg\n",action_H);
    
    // //compute action for G
    // const double action_G=
    //   gluonic_action(conf,theory_pars.gauge_action_name,theory_pars.beta);
    // VERBOSITY_LV1_MASTER_PRINTF("Gauge action: %lg\n",action_G);
    
    // return action_G+action_H;
    
    return {};
  }
  
  //perform a full hmc step and return the difference between final and original action
  double pure_gauge_hmc_step(quad_su3 *out_conf,quad_su3 *in_conf,theory_pars_t &theory_pars,pure_gauge_evol_pars_t &evol_pars,rat_approx_t *rat_exp_H,int itraj)
  {
    if(in_conf==out_conf) CRASH("in==out");
    
    //header
    double hmc_time=-take_time();
    MASTER_PRINTF("Trajectory %d (nmd: %d)\n",itraj,evol_pars.nmd_steps);
    
    //allocate the momenta and copy old conf into new
    quad_su3 *H=nissa_malloc("H",locVol,quad_su3);
    vector_copy(out_conf,in_conf);
    
    generate_hmc_momenta(H);
    
    //compute the action
    double init_action=pure_gauge_action(out_conf,theory_pars,evol_pars,H);
    VERBOSITY_LV2_MASTER_PRINTF("Init action: %lg\n",init_action);
    
    Omelyan_pure_gauge_evolver(H,out_conf,&theory_pars,&evol_pars);
    
    //compute the action
    double final_action=pure_gauge_action(out_conf,theory_pars,evol_pars,H);
    VERBOSITY_LV2_MASTER_PRINTF("Final action: %lg\n",final_action);
    
    //compute the diff
    double diff_action=final_action-init_action;
    
    //free stuff
    nissa_free(H);
    
    //take time
    hmc_time+=take_time();
    VERBOSITY_LV1_MASTER_PRINTF("Total time to perform pure gauge hmc step: %lg s\n",hmc_time);
    
    return diff_action;
  }
}
