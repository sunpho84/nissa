#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/momenta/MFACC.hpp"
#include "inverters/momenta/cg_invert_MFACC.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

//#define DEBUG

#include "hmc/momenta/momenta_action.hpp"
#include "operations/gauge_fixing.hpp"

namespace nissa
{
  //generate Fourier acceleration-related fields
  THREADABLE_FUNCTION_1ARG(generate_MFACC_fields, su3*,pi)
  {
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      herm_put_to_gauss(pi[ivol],&(loc_rnd_gen[ivol]),1);
    set_borders_invalid(pi);
  }
  THREADABLE_FUNCTION_END
  
  //compute the action for the Fourier acceleration-related fields - last bit of eq.4
  double MFACC_fields_action(su3 **phi)
  {
    //summ the square of pi
    double glb_action_lx[NDIM/2];
    for(int ifield=0;ifield<NDIM/2;ifield++)
      double_vector_glb_scalar_prod(&(glb_action_lx[ifield]),(double*)(phi[ifield]),(double*)(phi[ifield]),sizeof(su3)/sizeof(double)*loc_vol);
    
    return (glb_action_lx[0]+glb_action_lx[1])/2;
  }
  
  //Evolve Fourier acceleration related fields - eq.5
  THREADABLE_FUNCTION_5ARG(evolve_MFACC_fields, su3**,phi, quad_su3*,conf, double,kappa, su3**,pi, double,dt)
  {
    verbosity_lv2_master_printf("Evolving Fourier fields, dt=%lg\n",dt);
    
    //allocate
    su3 *F=nissa_malloc("temp",loc_vol+bord_vol,su3);
    su3 *temp=nissa_malloc("temp",loc_vol+bord_vol,su3);
    
    for(int ifield=0;ifield<NDIM/2;ifield++)
      {
#ifdef DEBUG
	double eps=1e-5;
	
	//store initial link and compute action
	su3 sto;
	su3_copy(sto,pi[ifield][0]);
	double act_ori=MFACC_momenta_action(pi,conf,kappa);
	
	//store derivative
	su3 nu_plus,nu_minus;
	su3_put_to_zero(nu_plus);
	su3_put_to_zero(nu_minus);
	
	for(int igen=0;igen<8;igen++)
	  {
	    //prepare increment and change
	    su3 ba;
	    su3_prod_double(ba,gell_mann_matr[igen],eps);
	    
	    //change -, compute action
	    su3_subt(pi[ifield][0],sto,ba);
	    double act_minus=MFACC_momenta_action(pi,conf,kappa);
	    
	    //change +, compute action
	    su3_summ(pi[ifield][0],sto,ba);
	    double act_plus=MFACC_momenta_action(pi,conf,kappa);
	    
	    //set back everything
	    su3_copy(phi[ifield][0],sto);
	    
	    printf("plus: %+016.016le, ori: %+016.016le, minus: %+016.016le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
	    double gr_plus=(act_plus-act_ori)/eps;
	    double gr_minus=(act_ori-act_minus)/eps;
	    su3_summ_the_prod_double(nu_plus,gell_mann_matr[igen],gr_plus/2);
	    su3_summ_the_prod_double(nu_minus,gell_mann_matr[igen],gr_minus/2);
	  }
	
	//take the average
	su3 nu;
	su3_summ(nu,nu_plus,nu_minus);
	su3_prodassign_double(nu,0.5);
#endif
	
        //compute
        apply_MFACC(temp,conf,kappa,pi[ifield]);
        apply_MFACC(F,conf,kappa,temp);
        
#ifdef DEBUG
	master_printf("Comparing MFACC fields derivative\n");
	master_printf("an\n");
	su3_print(F[0]);
	master_printf("nu\n");
	su3_print(nu);
	//crash("ciccio");
#endif
	//evolve
        double_vector_summ_double_vector_prod_double((double*)(phi[ifield]),(double*)(phi[ifield]),(double*)F,dt,loc_vol*sizeof(su3)/sizeof(double));
      }
    
    nissa_free(F);
    nissa_free(temp);
  }
  THREADABLE_FUNCTION_END
  
  //Evolve Fourier acceleration related momenta - eq.6
  THREADABLE_FUNCTION_3ARG(evolve_MFACC_momenta, su3**,pi, su3**,phi, double,dt)
  {
    verbosity_lv2_master_printf("Evolving Fourier momenta, dt=%lg\n",dt);
    
    for(int ifield=0;ifield<NDIM/2;ifield++)
      double_vector_summ_double_vector_prod_double((double*)(pi[ifield]),(double*)(pi[ifield]),(double*)(phi[ifield]),-dt,loc_vol*sizeof(su3)/sizeof(double));
    su3_print(pi[0][0]);
  }
  THREADABLE_FUNCTION_END
  
  //compute the QCD force originated from MFACC momenta (derivative of \pi^\dag MM \pi/2) w.r.t U
  THREADABLE_FUNCTION_5ARG(MFACC_momenta_QCD_force, quad_su3*,F, quad_su3*,conf, double,kappa, su3**,pi, bool,reset)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Computing QCD force originated by MFACC momenta (derivative of \\pi^\\dag MM \\pi/2) w.r.t U\n");
    
    su3 *temp=nissa_malloc("temp",loc_vol+bord_vol,su3);
    if(reset) vector_reset(F);
    
#ifdef DEBUG
    double eps=1e-5;
    
    //store initial link and compute action
    su3 sto;
    su3_copy(sto,conf[0][0]);
    double act_ori=MFACC_momenta_action(pi,conf,kappa);
    
    //store derivative
    su3 nu_plus,nu_minus;
    su3_put_to_zero(nu_plus);
    su3_put_to_zero(nu_minus);
    
    for(int igen=0;igen<8;igen++)
      {
	//prepare increment and change
	su3 ba;
	su3_prod_double(ba,gell_mann_matr[igen],eps/2);
	su3 exp_mod;
	safe_hermitian_exact_i_exponentiate(exp_mod,ba);
	
	//change -, compute action
	unsafe_su3_dag_prod_su3(conf[0][0],exp_mod,sto);
	double act_minus=MFACC_momenta_action(pi,conf,kappa);
	
	//change +, compute action
	unsafe_su3_prod_su3(conf[0][0],exp_mod,sto);
	double act_plus=MFACC_momenta_action(pi,conf,kappa);
	
	//set back everything
	su3_copy(conf[0][0],sto);
	
	//printf("plus: %+016.016le, ori: %+016.016le, minus: %+016.016le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
	double gr_plus=-(act_plus-act_ori)/eps;
	double gr_minus=-(act_ori-act_minus)/eps;
	su3_summ_the_prod_idouble(nu_plus,gell_mann_matr[igen],gr_plus/2);
	su3_summ_the_prod_idouble(nu_minus,gell_mann_matr[igen],gr_minus/2);
      }
    
    //take the average
    su3 nu;
    su3_summ(nu,nu_plus,nu_minus);
    su3_prodassign_double(nu,0.5);
    
    vector_reset(F);
#endif
    
    for(int ifield=0;ifield<NDIM/2;ifield++)
      {
	//apply
	apply_MFACC(temp,conf,kappa,pi[ifield]);
	
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      //temporary pieces
	      su3 t,E;
	      int up=loclx_neighup[ivol][mu];
	      
	      //forward piece
	      unsafe_su3_dag_prod_su3_dag(t,conf[ivol][mu],temp[ivol]);
	      unsafe_su3_prod_su3(E,pi[ifield][up],t);
	      
	      //backward piece
	      unsafe_su3_dag_prod_su3_dag(t,temp[up],conf[ivol][mu]);
	      su3_summ_the_prod_su3(E,t,pi[ifield][ivol]);
	      
	      //factor of kappa/16
	      su3_summ_the_prod_double(F[ivol][mu],E,-kappa/16);
	    }
	THREAD_BARRIER();
      }
    set_borders_invalid(F);
    
    nissa_free(temp);
    
#ifdef DEBUG
    su3 r1,r2;
    unsafe_su3_prod_su3(r1,conf[0][0],F[0][0]);
    unsafe_su3_traceless_anti_hermitian_part(r2,r1);
    
    master_printf("Comparing MFACC momenta QCD force\n");
    master_printf("an\n");
    su3_print(r2);
    master_printf("nu\n");
    su3_print(nu);
#endif
  }
  THREADABLE_FUNCTION_END
  
  //compute the QCD momenta force (derivative of \H^\dag G \H/2) (eq.8)
  THREADABLE_FUNCTION_7ARG(MFACC_QCD_momenta_QCD_force, quad_su3*,F, quad_su3*,conf, double,kappa, int,niter, double,residue, quad_su3*,H, bool,reset)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Computing QCD force due to Fourier Accelerated QCD momentas\n");
    
#ifdef DEBUG
    double eps=1e-5;
    
    //store initial link and compute action
    su3 sto;
    su3_copy(sto,conf[0][0]);
    double act_ori=momenta_action_with_FACC(conf,kappa,niter,residue,H);
    
    //store derivative
    su3 nu_plus,nu_minus;
    su3_put_to_zero(nu_plus);
    su3_put_to_zero(nu_minus);
    
    for(int igen=0;igen<8;igen++)
      {
	//prepare increment and change
	su3 ba;
	su3_prod_double(ba,gell_mann_matr[igen],eps/2);
	su3 exp_mod;
	safe_hermitian_exact_i_exponentiate(exp_mod,ba);
	
	//change -, compute action
	unsafe_su3_dag_prod_su3(conf[0][0],exp_mod,sto);
	double act_minus=momenta_action_with_FACC(conf,kappa,niter,residue,H);
	
	//change +, compute action
	unsafe_su3_prod_su3(conf[0][0],exp_mod,sto);
	double act_plus=momenta_action_with_FACC(conf,kappa,niter,residue,H);
	
	//set back everything
	su3_copy(conf[0][0],sto);
	
	//printf("plus: %+016.016le, ori: %+016.016le, minus: %+016.016le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
	double gr_plus=-(act_plus-act_ori)/eps;
	double gr_minus=-(act_ori-act_minus)/eps;
	su3_summ_the_prod_idouble(nu_plus,gell_mann_matr[igen],gr_plus/2);
	su3_summ_the_prod_idouble(nu_minus,gell_mann_matr[igen],gr_minus/2);
      }
    
    //take the average
    su3 nu;
    su3_summ(nu,nu_plus,nu_minus);
    su3_prodassign_double(nu,0.5);
    
    vector_reset(F);
#endif
    
    su3 *H_nu=nissa_malloc("H_nu",loc_vol+bord_vol,su3);
    su3 *temp=nissa_malloc("temp",loc_vol+bord_vol,su3);
    if(reset) vector_reset(F);
    
    for(int nu=0;nu<NDIM;nu++)
      {
	//copy out
        NISSA_PARALLEL_LOOP(ivol,0,loc_vol) su3_copy(H_nu[ivol],H[ivol][nu]);
        set_borders_invalid(H_nu);
	
	//invert
	inv_MFACC_cg(temp,NULL,conf,kappa,niter,residue,H_nu);
	
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  {
	    //D(A A^-1) = 0 = D(A) A^-1 + A D(A^-1) -> D(A^-1)=-A^-1 DA A^-1
	    
	    for(int mu=0;mu<NDIM;mu++)
	      {
		su3 E,t;
		int up=loclx_neighup[ivol][mu];
		
		//forward piece
		unsafe_su3_dag_prod_su3_dag(t,conf[ivol][mu],temp[ivol]);
		unsafe_su3_prod_su3(E,temp[up],t);
		
		//backward piece
		unsafe_su3_dag_prod_su3_dag(t,temp[up],conf[ivol][mu]);
		su3_summ_the_prod_su3(E,t,temp[ivol]);
		
		su3_summ_the_prod_double(F[ivol][mu],E,kappa/32); //it likes /32 more then -/32
	      }
	  }
	THREAD_BARRIER();
      }
    set_borders_invalid(F);
    
#ifdef DEBUG
    su3 r1,r2;
    unsafe_su3_prod_su3(r1,conf[0][0],F[0][0]);
    unsafe_su3_traceless_anti_hermitian_part(r2,r1);
    
    master_printf("QCD force originating from H\n");
    master_printf("an\n");
    su3_print(r2);
    master_printf("nu\n");
    su3_print(nu);
    //crash("ciccio");
#endif
    
    nissa_free(temp);
    nissa_free(H_nu);
  }
  THREADABLE_FUNCTION_END
}
