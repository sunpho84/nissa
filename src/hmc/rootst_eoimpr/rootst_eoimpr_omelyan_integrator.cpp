#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "hmc/backfield.hpp"
#include "new_types/complex.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "hmc/gauge/gluonic_force.hpp"
#include "hmc/gauge/topological_force.hpp"

#include "rootst_eoimpr_quark_force.hpp"

namespace nissa
{
  //unitarize the conf by explicitly inverting it
  THREADABLE_FUNCTION_1ARG(lx_conf_unitarize_maximal_trace_projecting, quad_su3*,conf)
  {
    GET_THREAD_ID();
    
    addrem_stagphases_to_lx_conf(conf);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<4;mu++)
	su3_unitarize_maximal_trace_projecting(conf[ivol][mu],conf[ivol][mu]);
    
    addrem_stagphases_to_lx_conf(conf);
    
    set_borders_invalid(conf);
  }
  THREADABLE_FUNCTION_END

  // Evolve momenta according to the pure gauge force
  // calculate H=H-F*dt to evolve link momenta
  // i.e calculate v(t+dt)=v(t)+a*dt
  THREADABLE_FUNCTION_5ARG(evolve_lx_momenta_with_pure_gauge_force, quad_su3*,H, quad_su3*,conf, theory_pars_t*,theory_pars, double,dt, quad_su3*,ext_F)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Evolving momenta with pure gauge force, dt=%lg\n",dt);
    
    //allocate force and compute it
    quad_su3 *F=(ext_F==NULL)?nissa_malloc("F",loc_vol,quad_su3):ext_F;
    compute_gluonic_force_lx_conf(F,conf,theory_pars);
    
    //evolve
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<4;mu++)
	for(int ic1=0;ic1<3;ic1++)
	  for(int ic2=0;ic2<3;ic2++)
	    complex_subt_the_prod_idouble(H[ivol][mu][ic1][ic2],F[ivol][mu][ic1][ic2],dt);
    
    if(ext_F==NULL) nissa_free(F);
    else THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
  
  //evolve the momenta with topological force
  THREADABLE_FUNCTION_5ARG(evolve_lx_momenta_with_topological_force, quad_su3*,H, quad_su3*,conf, topotential_pars_t*,topars, double,dt, quad_su3*,ext_F)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Evolving momenta with topological force, dt=%lg\n",dt);
    
    //allocate force and compute it
    quad_su3 *F=(ext_F==NULL)?nissa_malloc("F",loc_vol,quad_su3):ext_F;
    compute_topological_force_lx_conf(F,conf,topars);
    
    //evolve
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<4;mu++)
	for(int ic1=0;ic1<3;ic1++)
	  for(int ic2=0;ic2<3;ic2++)
	    complex_subt_the_prod_idouble(H[ivol][mu][ic1][ic2],F[ivol][mu][ic1][ic2],dt);
    
    if(ext_F==NULL) nissa_free(F);
    else THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END

  //evolve the configuration with the momenta
  THREADABLE_FUNCTION_4ARG(evolve_lx_conf_with_momenta, quad_su3*,lx_conf, quad_su3*,H, theory_pars_t*,tp, double,dt)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Evolving conf with momenta, dt=%lg\n",dt);
    
    //evolve
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<4;mu++)
	{
	  su3 t1,t2;
	  su3_prod_double(t1,H[ivol][mu],dt);
	  safe_anti_hermitian_exact_i_exponentiate(t2,t1);
	  
	  safe_su3_prod_su3(lx_conf[ivol][mu],t2,lx_conf[ivol][mu]);
	}
    
    set_borders_invalid(lx_conf);
  }
  THREADABLE_FUNCTION_END

  //evolve the configuration according to pure gauge
  THREADABLE_FUNCTION_4ARG(omelyan_pure_gauge_evolver_lx_conf, quad_su3*,H, quad_su3*,lx_conf, theory_pars_t*,theory_pars, hmc_evol_pars_t*,simul)
  {
    //macro step or micro step
    double dt=simul->traj_length/simul->nmd_steps/simul->ngauge_substeps/2,
      dth=dt/2,ldt=dt*OMELYAN_LAMBDA,l2dt=2*OMELYAN_LAMBDA*dt,uml2dt=(1-2*OMELYAN_LAMBDA)*dt;
    int nsteps=simul->ngauge_substeps;
    quad_su3 *aux_F=nissa_malloc("aux_F",loc_vol,quad_su3);
    
    topotential_pars_t *topars=&(theory_pars->topotential_pars);
    
    //     Compute H(t+lambda*dt) i.e. v1=v(t)+a[r(t)]*lambda*dt (first half step)
    evolve_lx_momenta_with_pure_gauge_force(H,lx_conf,theory_pars,ldt,aux_F);
    if(topars->flag) evolve_lx_momenta_with_topological_force(H,lx_conf,topars,ldt,aux_F);
    
    //         Main loop
    for(int istep=0;istep<nsteps;istep++)
      {
	verbosity_lv1_master_printf(" Omelyan gauge micro-step %d/%d\n",istep+1,nsteps);
	
	//decide if last step is final or not
	double last_dt=(istep==(nsteps-1)) ? ldt : l2dt;
	
	//     Compute U(t+dt/2) i.e. r1=r(t)+v1*dt/2
	evolve_lx_conf_with_momenta(lx_conf,H,theory_pars,dth);
	//     Compute H(t+(1-2*lambda)*dt) i.e. v2=v1+a[r1]*(1-2*lambda)*dt
	evolve_lx_momenta_with_pure_gauge_force(H,lx_conf,theory_pars,uml2dt,aux_F);
	if(topars->flag) evolve_lx_momenta_with_topological_force(H,lx_conf,topars,uml2dt,aux_F);
	//     Compute U(t+dt/2) i.e. r(t+dt)=r1+v2*dt/2
	evolve_lx_conf_with_momenta(lx_conf,H,theory_pars,dth);
	//     Compute H(t+dt) i.e. v(t+dt)=v2+a[r(t+dt)]*lambda*dt (at last step) or *2*lambda*dt
	evolve_lx_momenta_with_pure_gauge_force(H,lx_conf,theory_pars,last_dt,aux_F);
	if(topars->flag) evolve_lx_momenta_with_topological_force(H,lx_conf,topars,last_dt,aux_F);
	
	//normalize the configuration
	lx_conf_unitarize_maximal_trace_projecting(lx_conf);
      }
    
    nissa_free(aux_F);
  }
  THREADABLE_FUNCTION_END

  //wrapper
  void omelyan_pure_gauge_evolver_eo_conf(quad_su3 **H_eo,quad_su3 **conf_eo,theory_pars_t *theory_pars,hmc_evol_pars_t *simul)
  {
    quad_su3 *H_lx=nissa_malloc("H_lx",loc_vol,quad_su3);
    quad_su3 *conf_lx=nissa_malloc("conf_lx",loc_vol+bord_vol+edge_vol,quad_su3);
    
    paste_eo_parts_into_lx_conf(H_lx,H_eo);
    paste_eo_parts_into_lx_conf(conf_lx,conf_eo);
    
    omelyan_pure_gauge_evolver_lx_conf(H_lx,conf_lx,theory_pars,simul);
    
    split_lx_conf_into_eo_parts(H_eo,H_lx);
    split_lx_conf_into_eo_parts(conf_eo,conf_lx);
    
    nissa_free(conf_lx);
    nissa_free(H_lx);
  }

  /////////////////////////////////////// QUARK E/O PART ////////////////////////////////////////////////

  //unitarize the conf by explicitly inverting it
  THREADABLE_FUNCTION_1ARG(eo_conf_unitarize_maximal_trace_projecting, quad_su3**,conf)
  {
    GET_THREAD_ID();
    
    addrem_stagphases_to_eo_conf(conf);
    
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	  for(int mu=0;mu<4;mu++)
	    su3_unitarize_maximal_trace_projecting(conf[par][ivol][mu],conf[par][ivol][mu]);
	
	set_borders_invalid(conf[par]);
      }
    
    addrem_stagphases_to_eo_conf(conf);
  }
  THREADABLE_FUNCTION_END

  //evolve the b field
  void evolve_b_with_H_B_momenta(double *H_B,theory_pars_t *tp,double dt)
  {
    //evolve b
    if(H_B!=NULL) update_backfield(tp,tp->em_field_pars.B[tp->em_field_pars.meta.component]+(*H_B)*dt);
  }

  // Evolve momenta according to the rooted staggered force
  THREADABLE_FUNCTION_8ARG(evolve_momenta_with_quark_and_magnetic_force, quad_su3**,H, double*,H_B, quad_su3**,conf, color**,pf, theory_pars_t*,theory_pars, rat_approx_t*,appr, double,residue, double,dt)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Evolving momenta with quark force, dt=%lg\n",dt);
    
    //allocate forces
    quad_su3 *F[2]={nissa_malloc("F0",loc_volh,quad_su3),nissa_malloc("F1",loc_volh,quad_su3)};
    double *F_B=NULL;
    if(theory_pars->em_field_pars.flag==2) F_B=nissa_malloc("F_B",1,double);
    
    //compute the force
    compute_rootst_eoimpr_quark_and_magnetic_force(F,F_B,conf,pf,theory_pars,appr,residue);
    
    //evolve
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	  for(int mu=0;mu<4;mu++)
	    for(int ic1=0;ic1<3;ic1++)
	      for(int ic2=0;ic2<3;ic2++)
		complex_subt_the_prod_idouble(H[par][ivol][mu][ic1][ic2],F[par][ivol][mu][ic1][ic2],dt);
        nissa_free(F[par]);
      }
    
    //evolve the B momenta, if neeed
    if(H_B!=NULL && F_B!=NULL)
      {
	master_printf("Updating B mom from %lg, force %lg\n",*H_B,*F_B);
	if(IS_MASTER_THREAD) (*H_B)+=(*F_B)*dt;
	nissa_free(F_B);
      }
  }
  THREADABLE_FUNCTION_END

  ////////////////////////////////////// MACRO OMELYAN ////////////////////////////////////////////////

  THREADABLE_FUNCTION_7ARG(omelyan_rootst_eoimpr_evolver, quad_su3**,H, double*,H_B, quad_su3**,conf, color**,pf, theory_pars_t*,theory_pars, rat_approx_t*,appr, hmc_evol_pars_t*,simul)
  {
    //macro step or micro step
    double dt=simul->traj_length/simul->nmd_steps,
      ldt=dt*OMELYAN_LAMBDA,l2dt=2*OMELYAN_LAMBDA*dt,uml2dt=(1-2*OMELYAN_LAMBDA)*dt;  
    int nsteps=simul->nmd_steps;
    
    //     Compute H(t+lambda*dt) i.e. v1=v(t)+a[r(t)]*lambda*dt (first half step)
    evolve_momenta_with_quark_and_magnetic_force(H,H_B,conf,pf,theory_pars,appr,simul->md_residue,ldt);

    //         Main loop
    for(int istep=0;istep<nsteps;istep++)
      {
	verbosity_lv1_master_printf("Omelyan macro-step %d/%d\n",istep+1,nsteps);
	
	//decide if last step is final or not
	double last_dt=(istep==(nsteps-1)) ? ldt : l2dt;
	
	omelyan_pure_gauge_evolver_eo_conf(H,conf,theory_pars,simul);
	if(H_B) evolve_b_with_H_B_momenta(H_B,theory_pars,dt/2);
	evolve_momenta_with_quark_and_magnetic_force(H,H_B,conf,pf,theory_pars,appr,simul->md_residue,uml2dt);
	
	omelyan_pure_gauge_evolver_eo_conf(H,conf,theory_pars,simul);
	if(H_B) evolve_b_with_H_B_momenta(H_B,theory_pars,dt/2);
	evolve_momenta_with_quark_and_magnetic_force(H,H_B,conf,pf,theory_pars,appr,simul->md_residue,last_dt);
	
	//normalize the configuration
	eo_conf_unitarize_maximal_trace_projecting(conf);
      }
  }
  THREADABLE_FUNCTION_END
}
