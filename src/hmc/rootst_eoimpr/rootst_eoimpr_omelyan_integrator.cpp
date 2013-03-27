#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../geometry/geometry_eo.h"
#include "../../geometry/geometry_mix.h"
#include "../../new_types/complex.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../routines/ios.h"
#include "../../routines/openmp.h"

#include "../gauge/gluonic_force.h"

#include "rootst_eoimpr_quark_force.h"

//unitarize the conf by explicitly inverting it
THREADABLE_FUNCTION_2ARG(lx_conf_unitarize_explicitly_inverting, quad_su3*,conf, int,addrem_eo_phases)
{
  GET_THREAD_ID();
  
  if(addrem_eo_phases) crash("not implemented");
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<4;mu++)
      su3_unitarize_explicitly_inverting(conf[ivol][mu],conf[ivol][mu]);
  
  set_borders_invalid(conf);
}}

// Evolve momenta according to the pure gauge force
// calculate H=H-F*dt to evolve link momenta
// i.e calculate v(t+dt)=v(t)+a*dt
THREADABLE_FUNCTION_4ARG(evolve_lx_momenta_with_pure_gauge_force, quad_su3*,H, quad_su3*,conf, theory_pars_type*,theory_pars, double,dt)
{
  GET_THREAD_ID();
  
  verbosity_lv2_master_printf("Evolving momenta with force, dt=%lg\n",dt);
  
  //allocate force and compute it
  quad_su3 *F=nissa_malloc("F",loc_vol,quad_su3);
  compute_gluonic_force_lx_conf(F,conf,theory_pars);
  
  //evolve
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<4;mu++)
      for(int ic1=0;ic1<3;ic1++)
	for(int ic2=0;ic2<3;ic2++)
	  complex_subt_the_prod_idouble(H[ivol][mu][ic1][ic2],F[ivol][mu][ic1][ic2],dt);

  nissa_free(F);
}}

THREADABLE_FUNCTION_3ARG(evolve_lx_conf_with_momenta, quad_su3*,lx_conf, quad_su3*,H, double,dt)
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
}}

//evolve the configuration according to pure gauge
THREADABLE_FUNCTION_4ARG(omelyan_pure_gauge_evolver_lx_conf, quad_su3*,H, quad_su3*,lx_conf, theory_pars_type*,theory_pars, hmc_evol_pars_type*,simul)
{
  //macro step or micro step
  double dt=simul->traj_length/simul->nmd_steps/simul->ngauge_substeps/2,
    dth=dt/2,ldt=dt*OMELYAN_LAMBDA,l2dt=2*OMELYAN_LAMBDA*dt,uml2dt=(1-2*OMELYAN_LAMBDA)*dt;
  int nsteps=simul->ngauge_substeps;
  
  //     Compute H(t+lambda*dt) i.e. v1=v(t)+a[r(t)]*lambda*dt (first half step)
  evolve_lx_momenta_with_pure_gauge_force(H,lx_conf,theory_pars,ldt);
  
  //         Main loop
  for(int istep=0;istep<nsteps;istep++)
    {
      verbosity_lv1_master_printf(" Omelyan micro-step %d/%d\n",istep+1,nsteps);
      
      //decide if last step is final or not
      double last_dt=(istep==(nsteps-1)) ? ldt : l2dt;

      //     Compute U(t+dt/2) i.e. r1=r(t)+v1*dt/2
      evolve_lx_conf_with_momenta(lx_conf,H,dth);
      //     Compute H(t+(1-2*lambda)*dt) i.e. v2=v1+a[r1]*(1-2*lambda)*dt
      evolve_lx_momenta_with_pure_gauge_force(H,lx_conf,theory_pars,uml2dt);
      //     Compute U(t+dt/2) i.e. r(t+dt)=r1+v2*dt/2
      evolve_lx_conf_with_momenta(lx_conf,H,dth);
      //     Compute H(t+dt) i.e. v(t+dt)=v2+a[r(t+dt)]*lambda*dt (at last step) or *2*lambda*dt
      evolve_lx_momenta_with_pure_gauge_force(H,lx_conf,theory_pars,last_dt);
      
      //normalize the configuration
      lx_conf_unitarize_explicitly_inverting(lx_conf,false);
    }
}}

//wrapper
void omelyan_pure_gauge_evolver_eo_conf(quad_su3 **H_eo,quad_su3 **conf_eo,theory_pars_type *theory_pars,hmc_evol_pars_type *simul)
{
  quad_su3 *H_lx=nissa_malloc("H_lx",loc_vol,quad_su3);
  quad_su3 *conf_lx=nissa_malloc("conf_lx",loc_vol+bord_vol,quad_su3);
  
  paste_eo_parts_into_lx_conf(H_lx,H_eo);
  paste_eo_parts_into_lx_conf(conf_lx,conf_eo);
  
  omelyan_pure_gauge_evolver_lx_conf(H_lx,conf_lx,theory_pars,simul);
  
  nissa_free(conf_lx);
  nissa_free(H_lx);
}

/////////////////////////////////////// QUARK E/O PART ////////////////////////////////////////////////

//unitarize the conf by explicitly inverting it
THREADABLE_FUNCTION_2ARG(eo_conf_unitarize_explicitly_inverting, quad_su3**,conf, int,addrem_eo_phases)
{
  GET_THREAD_ID();
  
  if(addrem_eo_phases) addrem_stagphases_to_eo_conf(conf);
  
  for(int par=0;par<2;par++)
    {
      NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	for(int mu=0;mu<4;mu++)
	  su3_unitarize_explicitly_inverting(conf[par][ivol][mu],conf[par][ivol][mu]);
      
      set_borders_invalid(conf[par]);
    }
      
  if(addrem_eo_phases) addrem_stagphases_to_eo_conf(conf);
}}

//evolve the configuration by using the computed momenta
THREADABLE_FUNCTION_3ARG(evolve_eo_conf_with_momenta, quad_su3**,eo_conf, quad_su3**,H, double,dt)
{
  GET_THREAD_ID();
  
  verbosity_lv2_master_printf("Evolving conf with momenta, dt=%lg\n",dt);
  
  //evolve
  for(int par=0;par<2;par++)
    {
      NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	for(int mu=0;mu<4;mu++)
	  {
	    su3 t1,t2;
	    su3_prod_double(t1,H[par][ivol][mu],dt);
	    safe_anti_hermitian_exact_i_exponentiate(t2,t1);
	    
	    safe_su3_prod_su3(eo_conf[par][ivol][mu],t2,eo_conf[par][ivol][mu]);
	  }
      set_borders_invalid(eo_conf[par]);
    }
}}

// Evolve momenta according to the rooted staggered force
THREADABLE_FUNCTION_7ARG(evolve_momenta_with_quark_rootst_eoimpr_force, quad_su3**,H, quad_su3**,conf, color**,pf, theory_pars_type*,theory_pars, rat_approx_type*,appr, double,residue, double,dt)
{
  GET_THREAD_ID();
  
  verbosity_lv2_master_printf("Evolving momenta with force, dt=%lg\n",dt);
  
  //allocate force
  quad_su3 *F[2]={nissa_malloc("F0",loc_volh,quad_su3),nissa_malloc("F1",loc_volh,quad_su3)};
 
  //compute the force
  compute_rootst_eoimpr_quark_force(F,conf,pf,theory_pars,appr,residue);
 
  //evolve
  for(int par=0;par<2;par++)
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
      for(int mu=0;mu<4;mu++)
	for(int ic1=0;ic1<3;ic1++)
	  for(int ic2=0;ic2<3;ic2++)
	    complex_subt_the_prod_idouble(H[par][ivol][mu][ic1][ic2],F[par][ivol][mu][ic1][ic2],dt);

  for(int par=0;par<2;par++) nissa_free(F[par]);
}}

////////////////////////////////////// MACRO OMELYAN ////////////////////////////////////////////////

THREADABLE_FUNCTION_6ARG(omelyan_rootst_eoimpr_evolver, quad_su3**,H, quad_su3**,conf, color**,pf, theory_pars_type*,theory_pars, rat_approx_type*,appr, hmc_evol_pars_type*,simul)
{
  //macro step or micro step
  double dt=simul->traj_length/simul->nmd_steps,
    ldt=dt*OMELYAN_LAMBDA,l2dt=2*OMELYAN_LAMBDA*dt,uml2dt=(1-2*OMELYAN_LAMBDA)*dt;  
  int nsteps=simul->nmd_steps;
  
  //     Compute H(t+lambda*dt) i.e. v1=v(t)+a[r(t)]*lambda*dt (first half step)
  evolve_momenta_with_quark_rootst_eoimpr_force(H,conf,pf,theory_pars,appr,simul->md_residue,ldt);
  
  //         Main loop
  for(int istep=0;istep<nsteps;istep++)
    {
      verbosity_lv1_master_printf("Omelyan macro-step %d/%d\n",istep+1,nsteps);
      
      //decide if last step is final or not
      double last_dt=(istep==(nsteps-1)) ? ldt : l2dt;

      omelyan_pure_gauge_evolver_eo_conf(H,conf,theory_pars,simul);
      evolve_momenta_with_quark_rootst_eoimpr_force(H,conf,pf,theory_pars,appr,simul->md_residue,uml2dt);
      omelyan_pure_gauge_evolver_eo_conf(H,conf,theory_pars,simul);
      evolve_momenta_with_quark_rootst_eoimpr_force(H,conf,pf,theory_pars,appr,simul->md_residue,last_dt);
      
      //normalize the configuration
      eo_conf_unitarize_explicitly_inverting(conf,true);
    }
}}

