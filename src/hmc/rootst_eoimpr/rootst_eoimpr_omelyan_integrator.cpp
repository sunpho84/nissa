
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
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/smearing/stout.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "hmc/gauge/gluonic_force.hpp"
#include "hmc/gauge/topological_force.hpp"
#include "hmc/gauge/pure_gauge_omelyan_integrator.hpp"
#include "hmc/momenta/momenta_evolve.hpp"

#include "rootst_eoimpr_action.hpp"
#include "rootst_eoimpr_quark_force.hpp"

#define TOPO_MICRO 0
#define TOPO_MACRO 1

#define TOPO_EVOLUTION TOPO_MACRO

namespace nissa
{
  //evolve the momenta with topological force
  THREADABLE_FUNCTION_6ARG(evolve_lx_momenta_with_topological_force, quad_su3*,H, quad_su3*,conf, topotential_pars_t*,topars, double,dt, quad_su3*,ext_F, bool,phase_pres)
  {
    verbosity_lv2_master_printf("Evolving lx momenta with topological force, dt=%lg\n",dt);
    
    //allocate force and compute it
    quad_su3 *F=(ext_F==NULL)?nissa_malloc("F",loc_vol,quad_su3):ext_F;
    compute_topological_force_lx_conf(F,conf,topars,phase_pres);
    
    //evolve
    evolve_lx_momenta_with_force(H,F,dt);
    if(ext_F==NULL) nissa_free(F);
  }
  THREADABLE_FUNCTION_END
  
  //eo wrapper
  THREADABLE_FUNCTION_5ARG(evolve_eo_momenta_with_topological_force, quad_su3**,eo_H, quad_su3**,eo_conf, topotential_pars_t*,topars, double,dt, bool,phase_pres)
  {
    verbosity_lv2_master_printf("Evolving e/o momenta with topological force, dt=%lg\n",dt);
    
    //reorder
    quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+bord_vol+edge_vol,quad_su3);
    quad_su3 *lx_H=nissa_malloc("lx_H",loc_vol,quad_su3);
    paste_eo_parts_into_lx_vector(lx_conf,eo_conf);
    paste_eo_parts_into_lx_vector(lx_H,eo_H);
    
    evolve_lx_momenta_with_topological_force(lx_H,lx_conf,topars,dt,NULL,phase_pres);
    
    split_lx_vector_into_eo_parts(eo_H,lx_H);
    nissa_free(lx_H);
    nissa_free(lx_conf);
  }
  THREADABLE_FUNCTION_END
  
  //evolve the configuration according to pure gauge - note that there is a similar routine in "pure_gage"
  THREADABLE_FUNCTION_4ARG(omelyan_pure_gauge_evolver_lx_conf, quad_su3*,H, quad_su3*,lx_conf, theory_pars_t*,theory_pars, hmc_evol_pars_t*,simul)
  {
    addrem_stagphases_to_lx_conf(lx_conf);
    
    //macro step or micro step
    double dt=simul->traj_length/simul->nmd_steps/simul->ngauge_substeps/2,
      dth=dt/2,ldt=dt*OMELYAN_LAMBDA,l2dt=2*OMELYAN_LAMBDA*dt,uml2dt=(1-2*OMELYAN_LAMBDA)*dt;
    int nsteps=simul->ngauge_substeps;
    quad_su3 *aux_F=nissa_malloc("aux_F",loc_vol,quad_su3);
    
    topotential_pars_t *topars=&(theory_pars->topotential_pars);
    
    //     Compute H(t+lambda*dt) i.e. v1=v(t)+a[r(t)]*lambda*dt (first half step)
    evolve_momenta_with_pure_gauge_force(H,lx_conf,theory_pars,ldt,aux_F);
    if(topars->flag && TOPO_EVOLUTION==TOPO_MICRO) evolve_lx_momenta_with_topological_force(H,lx_conf,topars,ldt,aux_F,false);
    
    //         Main loop
    for(int istep=0;istep<nsteps;istep++)
      {
	verbosity_lv1_master_printf(" Omelyan gauge micro-step %d/%d\n",istep+1,nsteps);
	
	//decide if last step is final or not
	double last_dt=(istep==(nsteps-1)) ? ldt : l2dt;
	
	//     Compute U(t+dt/2) i.e. r1=r(t)+v1*dt/2
	evolve_lx_conf_with_momenta(lx_conf,H,dth);
	//     Compute H(t+(1-2*lambda)*dt) i.e. v2=v1+a[r1]*(1-2*lambda)*dt
	evolve_momenta_with_pure_gauge_force(H,lx_conf,theory_pars,uml2dt,aux_F);
	if(topars->flag && TOPO_EVOLUTION==TOPO_MICRO) evolve_lx_momenta_with_topological_force(H,lx_conf,topars,uml2dt,aux_F,false);
	//     Compute U(t+dt/2) i.e. r(t+dt)=r1+v2*dt/2
	evolve_lx_conf_with_momenta(lx_conf,H,dth);
	//     Compute H(t+dt) i.e. v(t+dt)=v2+a[r(t+dt)]*lambda*dt (at last step) or *2*lambda*dt
	evolve_momenta_with_pure_gauge_force(H,lx_conf,theory_pars,last_dt,aux_F);
	if(topars->flag && TOPO_EVOLUTION==TOPO_MICRO) evolve_lx_momenta_with_topological_force(H,lx_conf,topars,last_dt,aux_F,false);
	
	//normalize the configuration
	unitarize_lx_conf_maximal_trace_projecting(lx_conf);
      }
    
    addrem_stagphases_to_lx_conf(lx_conf);
    nissa_free(aux_F);
  }
  THREADABLE_FUNCTION_END
  
  //wrapper
  void omelyan_pure_gauge_evolver_eo_conf(quad_su3 **H_eo,quad_su3 **conf_eo,theory_pars_t *theory_pars,hmc_evol_pars_t *simul)
  {
    quad_su3 *H_lx=nissa_malloc("H_lx",loc_vol,quad_su3);
    quad_su3 *conf_lx=nissa_malloc("conf_lx",loc_vol+bord_vol+edge_vol,quad_su3);
    
    paste_eo_parts_into_lx_vector(H_lx,H_eo);
    paste_eo_parts_into_lx_vector(conf_lx,conf_eo);
    
    omelyan_pure_gauge_evolver_lx_conf(H_lx,conf_lx,theory_pars,simul);
    
    split_lx_vector_into_eo_parts(H_eo,H_lx);
    split_lx_vector_into_eo_parts(conf_eo,conf_lx);
    
    nissa_free(conf_lx);
    nissa_free(H_lx);
  }
  
  /////////////////////////////////////// QUARK E/O PART ////////////////////////////////////////////////
  
  // Evolve momenta according to the rooted staggered force
  THREADABLE_FUNCTION_6ARG(evolve_momenta_with_quark_force, quad_su3**,H, quad_su3**,conf, color***,pf, theory_pars_t*,theory_pars, hmc_evol_pars_t*,simul_pars, double,dt)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Evolving momenta with quark force, dt=%lg\n",dt);
    
    //allocate forces
    quad_su3 *F[2]={nissa_malloc("F0",loc_volh,quad_su3),nissa_malloc("F1",loc_volh,quad_su3)};
    
    //compute the force
    compute_rootst_eoimpr_quark_force(F,conf,pf,theory_pars,simul_pars->rat_appr,simul_pars->npseudo_fs,simul_pars->md_residue);
    
#if 0
    //print info on the norm of the force
    double norm2[2];
    for(int ieo=0;ieo<2;ieo++) double_vector_glb_scalar_prod(norm2+ieo,(double*)(F[ieo]),(double*)(F[ieo]),loc_volh*sizeof(quad_su3)/sizeof(double));
    master_printf("Fermionic force norm: %lg per site\n",sqrt(norm2[0]+norm2[1])/glb_vol);
    
    //pars for calc
    double eps=1e-4;
    int ivol=0,mu=0;
    
    //allocate smeared conf
    quad_su3 *sme_conf[2];
    for(int eo=0;eo<2;eo++) sme_conf[eo]=nissa_malloc("sme_conf",loc_volh+bord_volh+edge_volh,quad_su3);
    
    //smear
    addrem_stagphases_to_eo_conf(conf);
    stout_smear(sme_conf,conf,&(theory_pars->stout_pars));
    addrem_stagphases_to_eo_conf(conf);
    addrem_stagphases_to_eo_conf(sme_conf);
    
    //compute action before
    double act_bef;
    rootst_eoimpr_quark_action(&act_bef,sme_conf,theory_pars->nflavs,theory_pars->backfield,pf,simul_pars);
    
    //perform an inifinitesimal variation on site 0 dir 0
    su3 gen;
    su3_put_to_zero(gen);
    gen[1][0][0]=gen[0][1][0]=eps;
    su3 var;
    safe_anti_hermitian_exact_i_exponentiate(var,gen);
    
    su3_print(var);
    
    //modify
    su3 bef;
    su3_copy(bef,conf[0][ivol][mu]);
    if(rank==0 && IS_MASTER_THREAD) safe_su3_prod_su3(conf[0][ivol][mu],var,conf[0][ivol][mu]);
    set_borders_invalid(conf);
    
    //smear
    addrem_stagphases_to_eo_conf(conf);
    stout_smear(sme_conf,conf,&(theory_pars->stout_pars));
    addrem_stagphases_to_eo_conf(conf);
    addrem_stagphases_to_eo_conf(sme_conf);
    
    //compute action after rotation
    double act_aft;
    rootst_eoimpr_quark_action(&act_aft,sme_conf,theory_pars->nflavs,theory_pars->backfield,pf,simul_pars);
    
    //put back in place
    if(rank==0 && IS_MASTER_THREAD) su3_copy(conf[0][ivol][mu],bef);
    set_borders_invalid(conf);
    
    double f_num=(act_bef-act_aft)/eps;
    double f_ana=F[0][ivol][mu][1][0][IM]+F[0][ivol][mu][0][1][IM];
    master_printf("force: (%lg-%lg)/%lg=%lg numerical, %lg analytical\n",act_bef,act_aft,eps,f_num,f_ana);
    
    for(int eo=0;eo<2;eo++) nissa_free(sme_conf[eo]);
#endif
    
    //evolve
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	  for(int mu=0;mu<NDIM;mu++)
	    for(int ic1=0;ic1<NCOL;ic1++)
	      for(int ic2=0;ic2<NCOL;ic2++)
		complex_subt_the_prod_idouble(H[par][ivol][mu][ic1][ic2],F[par][ivol][mu][ic1][ic2],dt);
        nissa_free(F[par]);
      }
  }
  THREADABLE_FUNCTION_END
  
  ////////////////////////////////////// MACRO OMELYAN ////////////////////////////////////////////////
  
  THREADABLE_FUNCTION_5ARG(omelyan_rootst_eoimpr_evolver, quad_su3**,H, quad_su3**,conf, color***,pf, theory_pars_t*,theory_pars, hmc_evol_pars_t*,simul_pars)
  {
    //macro step or micro step
    double dt=simul_pars->traj_length/simul_pars->nmd_steps,
      ldt=dt*OMELYAN_LAMBDA,l2dt=2*OMELYAN_LAMBDA*dt,uml2dt=(1-2*OMELYAN_LAMBDA)*dt;
    int nsteps=simul_pars->nmd_steps;
    topotential_pars_t tp=theory_pars->topotential_pars;
    
    //     Compute H(t+lambda*dt) i.e. v1=v(t)+a[r(t)]*lambda*dt (first half step)
    evolve_momenta_with_quark_force(H,conf,pf,theory_pars,simul_pars,ldt);
    if(tp.flag && TOPO_EVOLUTION==TOPO_MACRO) evolve_eo_momenta_with_topological_force(H,conf,&tp,ldt,true);
    
    //         Main loop
    for(int istep=0;istep<nsteps;istep++)
      {
	verbosity_lv1_master_printf("Omelyan macro-step %d/%d\n",istep+1,nsteps);
	
	//decide if last step is final or not
	double last_dt=(istep==(nsteps-1)) ? ldt : l2dt;
	
	omelyan_pure_gauge_evolver_eo_conf(H,conf,theory_pars,simul_pars);
	evolve_momenta_with_quark_force(H,conf,pf,theory_pars,simul_pars,uml2dt);
	if(tp.flag && TOPO_EVOLUTION==TOPO_MACRO) evolve_eo_momenta_with_topological_force(H,conf,&tp,uml2dt,true);
	
	omelyan_pure_gauge_evolver_eo_conf(H,conf,theory_pars,simul_pars);
	evolve_momenta_with_quark_force(H,conf,pf,theory_pars,simul_pars,last_dt);
	if(tp.flag && TOPO_EVOLUTION==TOPO_MACRO) evolve_eo_momenta_with_topological_force(H,conf,&tp,last_dt,true);
	
	//normalize the configuration
	addrem_stagphases_to_eo_conf(conf);
	unitarize_eo_conf_maximal_trace_projecting(conf);
	addrem_stagphases_to_eo_conf(conf);
      }
  }
  THREADABLE_FUNCTION_END
}
