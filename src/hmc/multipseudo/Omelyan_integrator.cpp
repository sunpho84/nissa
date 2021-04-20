
#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <stdint.h>

#include "geometry/geometry_mix.hpp"
#include "operations/gaugeconf.hpp"

#include "hmc/hmc.hpp"
#include "hmc/theory_pars.hpp"
#include "hmc/gauge/gluonic_force.hpp"
#include "hmc/gauge/pure_gauge_Omelyan_integrator.hpp"
#include "hmc/gauge/topological_force.hpp"
#include "hmc/momenta/momenta_evolve.hpp"

#include "quark_force.hpp"
#include "theory_action.hpp"

#define TOPO_MICRO 0
#define TOPO_MACRO 1

#define TOPO_EVOLUTION TOPO_MACRO

namespace nissa
{
  //evolve the momenta with topological force
  void evolve_lx_momenta_with_topological_force(quad_su3* H,quad_su3* conf,topotential_pars_t* topars,double dt,quad_su3* ext_F)
  {
    verbosity_lv2_master_printf("Evolving lx momenta with topological force, dt=%lg\n",dt);
    
    //allocate force and compute it
    quad_su3 *F=(ext_F==NULL)?nissa_malloc("F",locVol.nastyConvert(),quad_su3):ext_F;
    compute_topological_force_lx_conf(F,conf,topars);
    
    //evolve
    evolve_lx_momenta_with_force(H,F,dt);
    if(ext_F==NULL) nissa_free(F);
  }
  
  //eo wrapper
  void evolve_eo_momenta_with_topological_force(eo_ptr<quad_su3> eo_H,eo_ptr<quad_su3> eo_conf,topotential_pars_t* topars,double dt)
  {
    verbosity_lv2_master_printf("Evolving e/o momenta with topological force, dt=%lg\n",dt);
    
    //reorder
    quad_su3 *lx_conf=nissa_malloc("lx_conf",locVolWithBordAndEdge.nastyConvert(),quad_su3);
    quad_su3 *lx_H=nissa_malloc("lx_H",locVol.nastyConvert(),quad_su3);
    paste_eo_parts_into_lx_vector(lx_conf,eo_conf);
    paste_eo_parts_into_lx_vector(lx_H,eo_H);
    
    evolve_lx_momenta_with_topological_force(lx_H,lx_conf,topars,dt,NULL);
    
    split_lx_vector_into_eo_parts(eo_H,lx_H);
    nissa_free(lx_H);
    nissa_free(lx_conf);
  }
  
  //evolve the configuration according to pure gauge - note that there is a similar routine in "pure_gage"
  void Omelyan_pure_gauge_evolver_lx_conf(quad_su3* H,quad_su3* lx_conf,theory_pars_t* theory_pars,hmc_evol_pars_t* simul)
  {
    //macro step or micro step
    double dt=simul->traj_length/simul->nmd_steps/simul->ngauge_substeps/2,
      dth=dt/2,ldt=dt*omelyan_lambda,l2dt=2*omelyan_lambda*dt,uml2dt=(1-2*omelyan_lambda)*dt;
    int nsteps=simul->ngauge_substeps;
    quad_su3 *aux_F=nissa_malloc("aux_F",locVol.nastyConvert(),quad_su3);
    
    topotential_pars_t *topars=&(theory_pars->topotential_pars);
    
    //     Compute H(t+lambda*dt) i.e. v1=v(t)+a[r(t)]*lambda*dt (first half step)
    evolve_momenta_with_pure_gauge_force(H,lx_conf,theory_pars,ldt,aux_F);
    if(topars->flag and TOPO_EVOLUTION==TOPO_MICRO) evolve_lx_momenta_with_topological_force(H,lx_conf,topars,ldt,aux_F);
    
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
	if(topars->flag and TOPO_EVOLUTION==TOPO_MICRO) evolve_lx_momenta_with_topological_force(H,lx_conf,topars,uml2dt,aux_F);
	//     Compute U(t+dt/2) i.e. r(t+dt)=r1+v2*dt/2
	evolve_lx_conf_with_momenta(lx_conf,H,dth);
	//     Compute H(t+dt) i.e. v(t+dt)=v2+a[r(t+dt)]*lambda*dt (at last step) or *2*lambda*dt
	evolve_momenta_with_pure_gauge_force(H,lx_conf,theory_pars,last_dt,aux_F);
	if(topars->flag and TOPO_EVOLUTION==TOPO_MICRO) evolve_lx_momenta_with_topological_force(H,lx_conf,topars,last_dt,aux_F);
	
	//normalize the configuration
	//unitarize_lx_conf_maximal_trace_projecting(lx_conf);
      }
    
    nissa_free(aux_F);
  }
  
  //wrapper
  void Omelyan_pure_gauge_evolver_eo_conf(eo_ptr<quad_su3> H_eo,eo_ptr<quad_su3> conf_eo,theory_pars_t *theory_pars,hmc_evol_pars_t *simul)
  {
    quad_su3 *H_lx=nissa_malloc("H_lx",locVol.nastyConvert(),quad_su3);
    quad_su3 *conf_lx=nissa_malloc("conf_lx",locVolWithBordAndEdge.nastyConvert(),quad_su3);
    
    paste_eo_parts_into_lx_vector(H_lx,H_eo);
    paste_eo_parts_into_lx_vector(conf_lx,conf_eo);
    
    Omelyan_pure_gauge_evolver_lx_conf(H_lx,conf_lx,theory_pars,simul);
    
    split_lx_vector_into_eo_parts(H_eo,H_lx);
    split_lx_vector_into_eo_parts(conf_eo,conf_lx);
    
    nissa_free(conf_lx);
    nissa_free(H_lx);
  }
  
  /////////////////////////////////////// QUARK E/O PART ////////////////////////////////////////////////
  
  // Evolve momenta according to the rooted staggered force
  void evolve_momenta_with_quark_force(eo_ptr<quad_su3> H,eo_ptr<quad_su3> conf,std::vector<std::vector<pseudofermion_t> >* pf,theory_pars_t* theory_pars,hmc_evol_pars_t* simul_pars,std::vector<rat_approx_t>* rat_appr,double dt)
  {
    
    verbosity_lv2_master_printf("Evolving momenta with quark force, dt=%lg\n",dt);
    
    //allocate forces
    eo_ptr<quad_su3> F={nissa_malloc("F0",locVolh.nastyConvert(),quad_su3),nissa_malloc("F1",locVolh.nastyConvert(),quad_su3)};
    
    //compute the force
    compute_quark_force(F,conf,pf,theory_pars,rat_appr,simul_pars->md_residue);
    
    //#define DEBUG_FORCE
    
#ifdef DEBUG_FORCE
    int par=0,ieo=1,mu=1;
    double eps=1e-4;
    
    //store initial link
    su3 sto;
    su3_copy(sto,conf[par][ieo][mu]);
    
    //allocate smeared conf
    std::array<quad_su3*,2> sme_conf;
    for(int eo=0;eo<2;eo++) sme_conf[eo]=nissa_malloc("sme_conf",loc_volh+bord_volh+edge_volh,quad_su3);
    
    //compute action before
    double act_ori;
    stout_smear(sme_conf,conf,&(theory_pars->stout_pars));
    compute_quark_action(&act_ori,sme_conf,theory_pars->backfield,pf,theory_pars->quarks,simul_pars,rat_appr);
    
    //store derivative
    su3 nu_plus,nu_minus;
    su3_put_to_zero(nu_plus);
    su3_put_to_zero(nu_minus);
    
    for(int igen=0;igen<NCOL*NCOL-1;igen++)
      {
	//prepare increment and change
	su3 ba;
	su3_prod_double(ba,gell_mann_matr[igen],eps/2);
	
	su3 exp_mod;
	safe_hermitian_exact_i_exponentiate(exp_mod,ba);
	
	//change -, compute action
	unsafe_su3_dag_prod_su3(conf[par][ieo][mu],exp_mod,sto);
	double act_minus;
	stout_smear(sme_conf,conf,&(theory_pars->stout_pars));
	compute_quark_action(&act_minus,sme_conf,theory_pars->backfield,pf,theory_pars->quarks,simul_pars,rat_appr);
	
	//change +, compute action
	unsafe_su3_prod_su3(conf[par][ieo][mu],exp_mod,sto);
	double act_plus;
	stout_smear(sme_conf,conf,&(theory_pars->stout_pars));
	compute_quark_action(&act_plus,sme_conf,theory_pars->backfield,pf,theory_pars->quarks,simul_pars,rat_appr);
	
	//set back everything
	su3_copy(conf[par][ieo][mu],sto);
	
	//printf("plus: %+016.016le, ori: %+016.016le, minus: %+016.016le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
	double gr_plus=-(act_plus-act_ori)/eps;
	double gr_minus=-(act_ori-act_minus)/eps;
	su3_summ_the_prod_idouble(nu_plus,gell_mann_matr[igen],gr_plus);
	su3_summ_the_prod_idouble(nu_minus,gell_mann_matr[igen],gr_minus);
      }
    
    //take the average
    su3 nu;
    su3_summ(nu,nu_plus,nu_minus);
    su3_prodassign_double(nu,0.5);
    
    master_printf("checking quark force\n");
    master_printf("an\n");
    su3_print(F[par][ieo][mu]);
    master_printf("nu\n");
    su3_print(nu);
    // master_printf("nu_plus\n");
    // su3_print(nu_plus);
    // master_printf("nu_minus\n");
    // su3_print(nu_minus);
    //crash("anna");
#endif
    
    //evolve
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  FOR_ALL_DIRS(mu)
	    for(int ic1=0;ic1<NCOL;ic1++)
	      for(int ic2=0;ic2<NCOL;ic2++)
		complex_subt_the_prod_idouble(H[par][ieo.nastyConvert()][mu.nastyConvert()][ic1][ic2],F[par][ieo.nastyConvert()][mu.nastyConvert()][ic1][ic2],dt);
	NISSA_PARALLEL_LOOP_END;
	
        nissa_free(F[par]);
      }
  }
  
  ////////////////////////////////////// MACRO OMELYAN ////////////////////////////////////////////////
  
  void Omelyan_integrator(eo_ptr<quad_su3> H,eo_ptr<quad_su3> conf,std::vector<std::vector<pseudofermion_t> >* pf,theory_pars_t* theory_pars,hmc_evol_pars_t* simul_pars,std::vector<rat_approx_t>* rat_appr)
  {
    int nsteps=simul_pars->nmd_steps;
    if(nsteps)
      {
	//macro step or micro step
	double dt=simul_pars->traj_length/simul_pars->nmd_steps,
	  ldt=dt*omelyan_lambda,l2dt=2*omelyan_lambda*dt,uml2dt=(1-2*omelyan_lambda)*dt;
	topotential_pars_t tp=theory_pars->topotential_pars;
	
	//     Compute H(t+lambda*dt) i.e. v1=v(t)+a[r(t)]*lambda*dt (first half step)
	evolve_momenta_with_quark_force(H,conf,pf,theory_pars,simul_pars,rat_appr,ldt);
	if(tp.flag and TOPO_EVOLUTION==TOPO_MACRO) evolve_eo_momenta_with_topological_force(H,conf,&tp,ldt);
	
	//         Main loop
	for(int istep=0;istep<nsteps;istep++)
	  {
	    verbosity_lv1_master_printf("Omelyan macro-step %d/%d\n",istep+1,nsteps);
	    
	    //decide if last step is final or not
	    double last_dt=(istep==(nsteps-1)) ? ldt : l2dt;
	    
	    Omelyan_pure_gauge_evolver_eo_conf(H,conf,theory_pars,simul_pars);
	    evolve_momenta_with_quark_force(H,conf,pf,theory_pars,simul_pars,rat_appr,uml2dt);
	    if(tp.flag and TOPO_EVOLUTION==TOPO_MACRO) evolve_eo_momenta_with_topological_force(H,conf,&tp,uml2dt);
	    
	    Omelyan_pure_gauge_evolver_eo_conf(H,conf,theory_pars,simul_pars);
	    evolve_momenta_with_quark_force(H,conf,pf,theory_pars,simul_pars,rat_appr,last_dt);
	    if(tp.flag and TOPO_EVOLUTION==TOPO_MACRO) evolve_eo_momenta_with_topological_force(H,conf,&tp,last_dt);
	  }
	
	//normalize the configuration
	unitarize_eo_conf_maximal_trace_projecting(conf);
      }
  }
}
