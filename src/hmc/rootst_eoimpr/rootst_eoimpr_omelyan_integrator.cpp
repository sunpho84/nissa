#include "../../base/global_variables.h"
#include "../../base/routines.h"
#include "../../base/vectors.h"
#include "../../geometry/geometry_eo.h"
#include "../../new_types/complex.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"

#include "../gauge/tree_level_Symanzik_action.h"

#include "rootst_eoimpr_force.h"


//unitarize the conf by explicitly inverting it
void eo_conf_unitarize_explicitly_inverting(quad_su3 **conf)
{
  addrem_stagphases_to_eo_conf(conf);
  
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  su3_unitarize_explicitly_inverting(conf[par][ivol][mu],conf[par][ivol][mu]);
      
      set_borders_invalid(conf[par]);
    }
      
  addrem_stagphases_to_eo_conf(conf);
}

// Evolve momenta according to the rooted staggered force
// calculate H=H-F*dt to evolve link momenta
// i.e calculate v(t+dt)=v(t)+a*dt
void evolve_momenta_with_full_rootst_eoimpr_force(quad_su3 **H,quad_su3 **conf,color **pf,theory_pars_type *theory_pars,rat_approx_type *appr,double residue,double dt,hmc_force_piece force_piece=BOTH_FORCE_PIECES)
{
  verbosity_lv2_master_printf("Evolving momenta with force, dt=%lg\n",dt);
  //allocate force
  quad_su3 *F[2]={nissa_malloc("F0",loc_volh,quad_su3),nissa_malloc("F1",loc_volh,quad_su3)};
  
  //compute the force
  full_rootst_eoimpr_force(F,conf,pf,theory_pars,appr,residue,force_piece);
  
  //evolve  
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  for(int ic1=0;ic1<3;ic1++)
	    for(int ic2=0;ic2<3;ic2++)
	      complex_subt_the_prod_idouble(H[par][ivol][mu][ic1][ic2],F[par][ivol][mu][ic1][ic2],dt);
      
      nissa_free(F[par]);
    }
}

//eolve the configuration by using the computed momenta
//this routine should be moved in a more general file
void evolve_conf_with_momenta(quad_su3 **eo_conf,quad_su3 **H,double dt)
{
  verbosity_lv2_master_printf("Evolving conf with momenta, dt=%lg\n",dt);
  
  //evolve
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  {
	    su3 t1,t2;
	    su3_prod_double(t1,H[par][ivol][mu],dt);
	    safe_anti_hermitian_exact_i_exponentiate(t2,t1);
	    
	    safe_su3_prod_su3(eo_conf[par][ivol][mu],t2,eo_conf[par][ivol][mu]);
	  }
      set_borders_invalid(eo_conf[par]);
    }
}

void omelyan_rootst_eoimpr_evolver(quad_su3 **H,quad_su3 **conf,color **pf,theory_pars_type *theory_pars,rat_approx_type *appr,hmc_evol_pars_type *simul,multistep_level multilev=MACRO_STEP)
{
  //define step length ant its multiples
  const double lambda=0.1931833;
  //macro step or micro step
  double dt;
  int nsteps;
  hmc_force_piece force_piece;
  if(simul->ngauge_substeps>1)
    {
      if(multilev==MACRO_STEP)
	{
	  dt=simul->traj_length/simul->nmd_steps;
	  force_piece=QUARK_FORCE_ONLY;
	  nsteps=simul->nmd_steps;
	}
      else
	{
	  dt=simul->traj_length/simul->nmd_steps/simul->ngauge_substeps/2;
	  force_piece=GAUGE_FORCE_ONLY;
	  nsteps=simul->ngauge_substeps;
	}
    }
  else
    {
      dt=simul->traj_length/simul->nmd_steps;
      force_piece=BOTH_FORCE_PIECES;
      nsteps=simul->nmd_steps;
    }
  double dth=dt/2;
  double ldt=dt*lambda,l2dt=2*lambda*dt,uml2dt=(1-2*lambda)*dt;  
  
  //     Compute H(t+lambda*dt) i.e. v1=v(t)+a[r(t)]*lambda*dt (first half step)
  evolve_momenta_with_full_rootst_eoimpr_force(H,conf,pf,theory_pars,appr,simul->md_residue,ldt,force_piece);
  
  //         Main loop
  for(int istep=0;istep<nsteps;istep++)
    {
      verbosity_lv1_master_printf("Omelyan step %d/%d\n",istep+1,nsteps);
      
      //decide if last step is final or not
      double last_dt=(istep==(nsteps-1)) ? ldt : l2dt;

      //     Compute U(t+dt/2) i.e. r1=r(t)+v1*dt/2
      if(multilev==MICRO_STEP||simul->ngauge_substeps<=1) evolve_conf_with_momenta(conf,H,dth);
      else omelyan_rootst_eoimpr_evolver(H,conf,pf,theory_pars,appr,simul,MICRO_STEP);
      //     Compute H(t+(1-2*lambda)*dt) i.e. v2=v1+a[r1]*(1-2*lambda)*dt
      evolve_momenta_with_full_rootst_eoimpr_force(H,conf,pf,theory_pars,appr,simul->md_residue,uml2dt,force_piece);
      //     Compute U(t+dt/2) i.e. r(t+dt)=r1+v2*dt/2
      if(multilev==MICRO_STEP||simul->ngauge_substeps<=1) evolve_conf_with_momenta(conf,H,dth);
      else omelyan_rootst_eoimpr_evolver(H,conf,pf,theory_pars,appr,simul,MICRO_STEP);
      //     Compute H(t+dt) i.e. v(t+dt)=v2+a[r(t+dt)]*lambda*dt (at last step) or *2*lambda*dt
      evolve_momenta_with_full_rootst_eoimpr_force(H,conf,pf,theory_pars,appr,simul->md_residue,last_dt,force_piece);
      
      //normalize the configuration
      eo_conf_unitarize_explicitly_inverting(conf);
    }
}
