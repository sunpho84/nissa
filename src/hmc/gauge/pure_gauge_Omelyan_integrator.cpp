#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/hmc.hpp"
#include "hmc/momenta/momenta_evolve.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "operations/gaugeconf.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

#include "gluonic_action.hpp"
#include "gluonic_force.hpp"

#include "pure_gauge_Omelyan_integrator.hpp"

//#define DEBUG

#ifdef DEBUG
namespace nissa
{
  double pure_gauge_action(quad_su3 *conf,theory_pars_t &theory_pars,pure_gauge_evol_pars_t &evol_pars,quad_su3 *H,su3 **phi,su3 **pi);
}
#endif

namespace nissa
{
  // Evolve momenta according to the pure gauge force
  // calculate H=H-F*dt to evolve link momenta
  // i.e calculate v(t+dt)=v(t)+a*dt
  void evolve_momenta_with_pure_gauge_force(LxField<quad_su3>& H,
					    const LxField<quad_su3>& conf,
					    const theory_pars_t& theory_pars,
					    const double& dt,
					    LxField<quad_su3>& F)
  {
    verbosity_lv2_master_printf("Evolving momenta with pure gauge force, dt=%lg\n",dt);
    
    compute_gluonic_force_lx_conf(F,conf,theory_pars);
    
    //evolve
    evolve_lx_momenta_with_force(H,F,dt);
  }
  
  //integrator for pure gauge
  void Omelyan_pure_gauge_evolver(quad_su3* H,quad_su3* conf,theory_pars_t* theory_pars,pure_gauge_evol_pars_t* simul)
  {
    //macro step or micro step
    double dt=simul->traj_length/simul->nmd_steps,dth=dt/2,ldt=dt*omelyan_lambda,l2dt=2*omelyan_lambda*dt,uml2dt=(1-2*omelyan_lambda)*dt;
    int nsteps=simul->nmd_steps;
    
    quad_su3 *F=nissa_malloc("F",locVol,quad_su3);
    
    //first evolve for momenta
    evolve_momenta_with_pure_gauge_force(H,conf,theory_pars,ldt,F);
    
    //         Main loop
    for(int istep=0;istep<nsteps;istep++)
      {
	verbosity_lv1_master_printf("Omelyan pure gauge step %d/%d\n",istep+1,nsteps);
	
	//decide if last step is final or not
	double last_dt=(istep==(nsteps-1)) ? ldt : l2dt;
	
	evolve_lx_conf_with_momenta(conf,H,dth);
	evolve_momenta_with_pure_gauge_force(H,conf,theory_pars,uml2dt,F);
	
	evolve_lx_conf_with_momenta(conf,H,dth);
	evolve_momenta_with_pure_gauge_force(H,conf,theory_pars,last_dt,F);
	
	//normalize the configuration
	unitarize_lx_conf_maximal_trace_projecting(conf);
      }
    
    nissa_free(F);
  }
}
