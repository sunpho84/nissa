#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/hmc.hpp"
#include "hmc/momenta/momenta_evolve.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "operations/gaugeconf.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "gluonic_action.hpp"
#include "gluonic_force.hpp"
#include "MFACC_fields.hpp"

#include "pure_gauge_Omelyan_integrator.hpp"

int evolve_FACC=1;

namespace nissa
{
  // Evolve momenta according to the pure gauge force
  // calculate H=H-F*dt to evolve link momenta
  // i.e calculate v(t+dt)=v(t)+a*dt
  THREADABLE_FUNCTION_5ARG(evolve_momenta_with_pure_gauge_force, quad_su3*,H, quad_su3*,conf, theory_pars_t*,theory_pars, double,dt, quad_su3*,ext_F)
  {
    verbosity_lv2_master_printf("Evolving momenta with pure gauge force, dt=%lg\n",dt);
    
    //allocate force and compute it
    quad_su3 *F=(ext_F==NULL)?nissa_malloc("F",loc_vol,quad_su3):ext_F;
    compute_gluonic_force_lx_conf(F,conf,theory_pars);
    
    //evolve
    evolve_lx_momenta_with_force(H,F,dt);
    
    if(ext_F==NULL) nissa_free(F);
  }
  THREADABLE_FUNCTION_END
  
  //same but with acceleration
  THREADABLE_FUNCTION_8ARG(evolve_momenta_and_FACC_momenta, quad_su3*,H, su3**,pi, quad_su3*,conf, su3**,phi, theory_pars_t*,theory_pars, pure_gauge_evol_pars_t*,simul, double,dt, quad_su3*,ext_F)
  {
    verbosity_lv2_master_printf("Evolving momenta and FACC momenta, dt=%lg\n",dt);
    
    quad_su3 *F=(ext_F==NULL)?nissa_malloc("F",loc_vol,quad_su3):ext_F;
    
    //evolve FACC momenta
    if(evolve_FACC&1) evolve_MFACC_momenta(pi,phi,dt);
    
    //compute the various contribution to the QCD force
    vector_reset(F);
    if(evolve_FACC&2) compute_gluonic_force_lx_conf(F,conf,theory_pars);
    if(evolve_FACC&1) summ_the_MFACC_momenta_QCD_force(F,conf,simul->kappa,pi);
    if(evolve_FACC&2) summ_the_MFACC_QCD_momenta_QCD_force(F,conf,simul->kappa,100000,simul->residue,H);
    evolve_lx_momenta_with_force(H,F,dt);
    
    if(ext_F==NULL) nissa_free(F);
  }
  THREADABLE_FUNCTION_END
  
  //combine the two fields evolution
  void evolve_lx_conf_with_accelerated_momenta_and_FACC_fields(quad_su3 *conf,su3 **phi,quad_su3 *H,su3 **pi,double kappa,int niter,double residue,double dt)
  {
    if(evolve_FACC&1) evolve_MFACC_fields(phi,conf,kappa,pi,dt);
    if(evolve_FACC&2) evolve_lx_conf_with_accelerated_momenta(conf,H,kappa,niter,residue,dt);
  }
  
  //integrator for pure gauge
  THREADABLE_FUNCTION_6ARG(Omelyan_pure_gauge_FACC_evolver, quad_su3*,H, quad_su3*,conf, su3**,pi, su3**,phi, theory_pars_t*,theory_pars, pure_gauge_evol_pars_t*,simul)
  {
    const int niter_max=100000;
    
    //macro step or micro step
    double dt=simul->traj_length/simul->nmd_steps,dth=dt/2,
      ldt=dt*OMELYAN_LAMBDA,l2dt=2*OMELYAN_LAMBDA*dt,uml2dt=(1-2*OMELYAN_LAMBDA)*dt;
    int nsteps=simul->nmd_steps;
    
    quad_su3 *F=nissa_malloc("F",loc_vol,quad_su3);
    
    evolve_momenta_and_FACC_momenta(H,pi,conf,phi,theory_pars,simul,ldt,F);
    
    //         Main loop
    for(int istep=0;istep<nsteps;istep++)
      {
	verbosity_lv1_master_printf("Omelyan step %d/%d\n",istep+1,nsteps);
	
	//decide if last step is final or not
	double last_dt=(istep==(nsteps-1)) ? ldt : l2dt;
	
	evolve_lx_conf_with_accelerated_momenta_and_FACC_fields(conf,phi,H,pi,simul->kappa,niter_max,simul->residue,dth);
	evolve_momenta_and_FACC_momenta(H,pi,conf,phi,theory_pars,simul,uml2dt,F);
	
	evolve_lx_conf_with_accelerated_momenta_and_FACC_fields(conf,phi,H,pi,simul->kappa,niter_max,simul->residue,dth);
	evolve_momenta_and_FACC_momenta(H,pi,conf,phi,theory_pars,simul,last_dt,F);
	
	//normalize the configuration
	unitarize_lx_conf_maximal_trace_projecting(conf);
      }
    
    nissa_free(F);
  }
  THREADABLE_FUNCTION_END
  
  //integrator for pure gauge
  THREADABLE_FUNCTION_4ARG(Omelyan_pure_gauge_evolver, quad_su3*,H, quad_su3*,conf, theory_pars_t*,theory_pars, pure_gauge_evol_pars_t*,simul)
  {
    //macro step or micro step
    double dt=simul->traj_length/simul->nmd_steps,dth=dt/2,ldt=dt*OMELYAN_LAMBDA,l2dt=2*OMELYAN_LAMBDA*dt,uml2dt=(1-2*OMELYAN_LAMBDA)*dt;
    int nsteps=simul->nmd_steps;
    
    quad_su3 *F=nissa_malloc("F",loc_vol,quad_su3);
    
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
  THREADABLE_FUNCTION_END
}
