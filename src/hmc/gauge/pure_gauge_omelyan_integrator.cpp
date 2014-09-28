
#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/momenta/momenta_evolve.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "operations/gaugeconf.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "gluonic_action.hpp"
#include "gluonic_force.hpp"
#include "MFACC_fields.hpp"

namespace nissa
{
  // Evolve momenta according to the pure gauge force
  // calculate H=H-F*dt to evolve link momenta
  // i.e calculate v(t+dt)=v(t)+a*dt
  THREADABLE_FUNCTION_6ARG(evolve_momenta_with_pure_gauge_force, quad_su3*,H, quad_su3*,conf, theory_pars_t*,theory_pars, double,dt, bool,phase_pres, quad_su3*,ext_F)
  {
    verbosity_lv2_master_printf("Evolving momenta with pure gauge force, dt=%lg\n",dt);
    
    //allocate force and compute it
    quad_su3 *F=(ext_F==NULL)?nissa_malloc("F",loc_vol,quad_su3):ext_F;
    compute_gluonic_force_lx_conf(F,conf,theory_pars,phase_pres);
    
    //evolve
    evolve_lx_momenta_with_force(H,F,dt);
    if(ext_F==NULL) nissa_free(F);
  }
  THREADABLE_FUNCTION_END
  
  //integrator for pure gauge - WARNING: assuming no phase present!
  THREADABLE_FUNCTION_6ARG(omelyan_pure_gauge_evolver, quad_su3*,H, quad_su3*,conf, su3**,pi, su3**,phi, theory_pars_t*,theory_pars, pure_gauge_evol_pars_t*,simul)
  {
    //macro step or micro step
    double dt=simul->traj_length/simul->nmd_steps,dth=dt/2,
      ldt=dt*OMELYAN_LAMBDA,l2dt=2*OMELYAN_LAMBDA*dt,uml2dt=(1-2*OMELYAN_LAMBDA)*dt;  
    int nsteps=simul->nmd_steps;
    
    double kappa=0.9;
    int niter=10000;
    double residue=1.e-10;
    
    //     Compute H(t+lambda*dt) i.e. v1=v(t)+a[r(t)]*lambda*dt (first half step)
    //evolve_momenta_with_pure_gauge_force(H,conf,theory_pars,ldt,false/*assuming no phase is present*/,NULL);
    if(pi!=NULL) evolve_MFACC_momenta(pi,phi,ldt);
    
    //         Main loop
    for(int istep=0;istep<nsteps;istep++)
      {
	verbosity_lv1_master_printf("Omelyan step %d/%d\n",istep+1,nsteps);
	
	//decide if last step is final or not
	double last_dt=(istep==(nsteps-1)) ? ldt : l2dt;
	
	//evolve_lx_conf_with_accelerated_momenta(conf,H,kappa,niter,residue,dth);
	if(phi!=NULL) evolve_MFACC_fields(phi,conf,kappa,pi,dth);
	//evolve_momenta_with_pure_gauge_force(H,conf,theory_pars,uml2dt,false,NULL);
	if(pi!=NULL) evolve_MFACC_momenta(pi,phi,uml2dt);
	
	//evolve_lx_conf_with_accelerated_momenta(conf,H,kappa,niter,residue,dth);
	if(phi!=NULL) evolve_MFACC_fields(phi,conf,kappa,pi,dth);
	//evolve_momenta_with_pure_gauge_force(H,conf,theory_pars,last_dt,false,NULL);
	if(pi!=NULL) evolve_MFACC_momenta(pi,phi,last_dt);
	
	//normalize the configuration
	unitarize_lx_conf_maximal_trace_projecting(conf);
      }
  }
  THREADABLE_FUNCTION_END
}
