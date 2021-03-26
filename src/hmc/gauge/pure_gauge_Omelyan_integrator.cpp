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
#include "MFACC_fields.hpp"

#include "pure_gauge_Omelyan_integrator.hpp"

//#define DEBUG

#ifdef DEBUG
namespace nissa
{
  double pure_gauge_action(quad_su3 *conf,theory_pars_t &theory_pars,pure_gauge_evol_pars_t &evol_pars,quad_su3 *H,su3 **phi,su3 **pi);
}
#endif

int evolve_SU3=1;
int evolve_FACC=1;

namespace nissa
{
  // Evolve momenta according to the pure gauge force
  // calculate H=H-F*dt to evolve link momenta
  // i.e calculate v(t+dt)=v(t)+a*dt
  void evolve_momenta_with_pure_gauge_force(quad_su3* H,quad_su3* conf,theory_pars_t* theory_pars,double dt,quad_su3* ext_F)
  {
    verbosity_lv2_master_printf("Evolving momenta with pure gauge force, dt=%lg\n",dt);
    
    //allocate force and compute it
    quad_su3 *F=(ext_F==NULL)?nissa_malloc("F",locVol,quad_su3):ext_F;
    compute_gluonic_force_lx_conf(F,conf,theory_pars);
    
    //evolve
    evolve_lx_momenta_with_force(H,F,dt);
    
    if(ext_F==NULL) nissa_free(F);
  }
  
  //same but with acceleration
  void evolve_momenta_and_FACC_momenta(quad_su3* H,su3** pi,quad_su3* conf,su3** phi,theory_pars_t* theory_pars,pure_gauge_evol_pars_t* simul,double dt,quad_su3* ext_F)
  {
    verbosity_lv2_master_printf("Evolving momenta and FACC momenta, dt=%lg\n",dt);
    
    quad_su3 *F=(ext_F==NULL)?nissa_malloc("F",locVol,quad_su3):ext_F;
    
#ifdef DEBUG
    vector_reset(F);
    double eps=1e-5;
    
    //store initial link and compute action
    su3 sto;
    su3_copy(sto,conf[0][0]);
    double act_ori=pure_gauge_action(conf,*theory_pars,*simul,H,phi,pi);
    
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
	unsafe_su3_dag_prod_su3(conf[0][0],exp_mod,sto);
	double act_minus=pure_gauge_action(conf,*theory_pars,*simul,H,phi,pi);
	
	//change +, compute action
	unsafe_su3_prod_su3(conf[0][0],exp_mod,sto);
	double act_plus=pure_gauge_action(conf,*theory_pars,*simul,H,phi,pi);
	
	//set back everything
	su3_copy(conf[0][0],sto);
	
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
    
    vector_reset(F);
#endif
    
    //compute the various contribution to the QCD force
    if(evolve_SU3)
      {
	//compute without TA
	vector_reset(F);
	compute_gluonic_force_lx_conf_do_not_finish(F,conf,theory_pars);
	summ_the_MFACC_momenta_QCD_force(F,conf,simul->kappa,pi,simul->naux_fields);
	summ_the_MFACC_QCD_momenta_QCD_force(F,conf,simul->kappa,100000,simul->residue,H);
	
	//finish the calculation
	gluonic_force_finish_computation(F,conf);
	
	evolve_lx_momenta_with_force(H,F,dt);
      }
    
#ifdef DEBUG
    master_printf("checking TOTAL gauge force\n");
    master_printf("an\n");
    su3_print(F[0][0]);
    master_printf("nu\n");
    su3_print(nu);
    master_printf("nu_plus\n");
    su3_print(nu_plus);
    master_printf("nu_minus\n");
    su3_print(nu_minus);
    //crash("anna");
#endif
    
    //evolve FACC momenta
    if(evolve_FACC) evolve_MFACC_momenta(pi,phi,simul->naux_fields,dt);
    
    if(ext_F==NULL) nissa_free(F);
  }
  
  //combine the two fields evolution
  void evolve_lx_conf_with_accelerated_momenta_and_FACC_fields(quad_su3 *conf,su3 **phi,quad_su3 *H,su3 **pi,int naux_fields,double kappa,int niter,double residue,double dt)
  {
    if(evolve_FACC) evolve_MFACC_fields(phi,naux_fields,conf,kappa,pi,dt);
    if(evolve_SU3) evolve_lx_conf_with_accelerated_momenta(conf,conf,H,kappa,niter,residue,dt);
  }
  
  //integrator for pure gauge
  void Omelyan_pure_gauge_FACC_evolver(quad_su3* H,quad_su3* conf,su3** pi,su3** phi,theory_pars_t* theory_pars,pure_gauge_evol_pars_t* simul)
  {
    const int niter_max=1000000;
    
    //macro step or micro step
    double dt=simul->traj_length/simul->nmd_steps,dth=dt/2,
      ldt=dt*omelyan_lambda,l2dt=2*omelyan_lambda*dt,uml2dt=(1-2*omelyan_lambda)*dt;
    int nsteps=simul->nmd_steps;
    
    quad_su3 *F=nissa_malloc("F",locVol,quad_su3);
    
    evolve_momenta_and_FACC_momenta(H,pi,conf,phi,theory_pars,simul,ldt,F);
    
    //         Main loop
    for(int istep=0;istep<nsteps;istep++)
      {
	verbosity_lv1_master_printf("Omelyan step %d/%d\n",istep+1,nsteps);
	
	//decide if last step is final or not
	double last_dt=(istep==(nsteps-1)) ? ldt : l2dt;
	
	evolve_lx_conf_with_accelerated_momenta_and_FACC_fields(conf,phi,H,pi,simul->naux_fields,simul->kappa,niter_max,simul->residue,dth);
	evolve_momenta_and_FACC_momenta(H,pi,conf,phi,theory_pars,simul,uml2dt,F);
	
	evolve_lx_conf_with_accelerated_momenta_and_FACC_fields(conf,phi,H,pi,simul->naux_fields,simul->kappa,niter_max,simul->residue,dth);
	evolve_momenta_and_FACC_momenta(H,pi,conf,phi,theory_pars,simul,last_dt,F);
	
	//normalize the configuration
	unitarize_lx_conf_maximal_trace_projecting(conf);
      }
    
    nissa_free(F);
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
