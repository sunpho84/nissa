#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/random.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/momenta/MFACC.hpp"
#include "inverters/momenta/cg_invert_MFACC.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#define DEBUG

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
	
	su3 nu_plus,nu_minus;
	for(int ic1=0;ic1<NCOL;ic1++)
	  for(int ic2=0;ic2<NCOL;ic2++)
	    for(int ri=0;ri<2;ri++)
	      {
		//prepare increment and change
		su3 ba;
		su3_put_to_zero(ba);
		ba[ic1][ic2][ri]+=eps/2;
		ba[ic2][ic1][ri]+=(ri==0?(+1):(-1))*eps/2;
		
		//change -, compute action
		su3_subt(pi[ifield][0],sto,ba);
		double act_minus=MFACC_momenta_action(pi,conf,kappa);
		
		//change +, compute action
		su3_summ(pi[ifield][0],sto,ba);
		double act_plus=MFACC_momenta_action(pi,conf,kappa);
		
		//set back everything
		su3_copy(phi[ifield][0],sto);
		
		//print pre and post, and write numerical
		//printf("plus: %+016.016le, ori: %+016.016le, minus: %+016.016le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
		nu_plus[ic1][ic2][ri]=-(act_plus-act_ori)/eps;
		nu_minus[ic1][ic2][ri]=-(act_ori-act_minus)/eps;
	      }
	
	//take the average
	su3 nu;
	su3_summ(nu,nu_plus,nu_minus);
#endif
	
        //compute
        apply_MFACC(temp,conf,kappa,pi[ifield]);
        apply_MFACC(F,conf,kappa,temp);
        
#ifdef DEBUG
	master_printf("an\n");
	su3_print(F[0]);
	master_printf("nu\n");
	su3_print(nu);
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
    
    for(int ifield=0;ifield<2;ifield++)
      double_vector_summ_double_vector_prod_double((double*)(pi[ifield]),(double*)(pi[ifield]),(double*)(phi[ifield]),-dt,loc_vol*sizeof(su3)/sizeof(double));
  }
  THREADABLE_FUNCTION_END
  
  //compute the MFACC momenta-related QCD force (derivative of \pi^\dag MM \pi/2)
  THREADABLE_FUNCTION_5ARG(MFACC_momenta_QCD_force, quad_su3*,F, quad_su3*,conf, double,kappa, su3**,pi, bool,reset)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Computing Fourier QCD force originated by acceleration-momenta\n");
    
    su3 *temp=nissa_malloc("temp",loc_vol+bord_vol,su3);
    if(reset) vector_reset(F);
    
#ifdef DEBUG
    vector_reset(F);
    double eps=1e-5;
    
    //store initial link and compute action
    su3 sto;
    su3_copy(sto,conf[0][0]);
    double act_ori=MFACC_momenta_action(pi,conf,kappa);
    
    //prepare increment and change
    su3 ba;
    su3_put_to_zero(ba);
    ba[1][0][0]=ba[0][1][0]=eps/2;
    su3 exp_mod;
    safe_anti_hermitian_exact_i_exponentiate(exp_mod,ba);
    su3_print(exp_mod);
    
    //change -, compute action
    unsafe_su3_dag_prod_su3(conf[0][0],exp_mod,sto);
    double act_minus=MFACC_momenta_action(pi,conf,kappa);
    
    //change +, compute action
    unsafe_su3_prod_su3(conf[0][0],exp_mod,sto);
    double act_plus=MFACC_momenta_action(pi,conf,kappa);
    
    //set back everything
    su3_copy(conf[0][0],sto);
    
    //print pre and post, and write numerical
    printf("plus: %+016.016le, ori: %+016.016le, minus: %+016.016le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
    double nu_plus=-(act_plus-act_ori)/eps;
    double nu_minus=-(act_ori-act_minus)/eps;
    double nu=-(act_plus-act_minus)/(2*eps);
    
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
	      su3_summ_the_prod_double(F[ivol][mu],E,kappa/16);
	    }
	THREAD_BARRIER();
      }
    set_borders_invalid(F);
    
    nissa_free(temp);
    
#ifdef DEBUG
    su3 r1,r2;
    unsafe_su3_prod_su3(r1,conf[0][0],F[0][0]);
    unsafe_su3_traceless_anti_hermitian_part(r2,r1);
    double tr=(r2[1][0][1]+r2[0][1][1])/2;
    printf("an: %+016.016le, nu: %+016.016le, nu+: %+016.016le, nu-: %+016.016le\n",tr,nu,nu_plus,nu_minus);
    //crash("anna");
#endif
  }
  THREADABLE_FUNCTION_END
  
  //compute the MFACC fields-related QCD force (derivative of \H^\dag G \H/2) (eq.8)
  THREADABLE_FUNCTION_7ARG(MFACC_QCD_momenta_QCD_force, quad_su3*,F, quad_su3*,conf, double,kappa, int,niter, double,residue, quad_su3*,H, bool,reset)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Computing Fourier acceleration fields originated QCD force\n");
    
    su3 *H_nu=nissa_malloc("H_nu",loc_vol+bord_vol,su3);
    su3 *temp=nissa_malloc("temp",loc_vol+bord_vol,su3);
    if(reset) vector_reset(F);
    
#ifdef DEBUG
    vector_reset(F);
    double eps=1e-5;
    
    //store initial link and compute action
    su3 sto;
    su3_copy(sto,conf[0][0]);
    double act_ori=momenta_action_with_FACC(conf,kappa,niter,residue,H);
    
    //prepare increment and change
    su3 ba;
    su3_put_to_zero(ba);
    ba[1][0][0]=ba[0][1][0]=eps/2;
    su3 exp_mod;
    safe_anti_hermitian_exact_i_exponentiate(exp_mod,ba);
    su3_print(exp_mod);
    
    //change -, compute action
    unsafe_su3_dag_prod_su3(conf[0][0],exp_mod,sto);
    double act_minus=momenta_action_with_FACC(conf,kappa,niter,residue,H);
    
    //change +, compute action
    unsafe_su3_prod_su3(conf[0][0],exp_mod,sto);
    double act_plus=momenta_action_with_FACC(conf,kappa,niter,residue,H);
    
    //set back everything
    su3_copy(conf[0][0],sto);
    
    //print pre and post, and write numerical
    printf("plus: %+016.016le, ori: %+016.016le, minus: %+016.016le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
    double nu_plus=-(act_plus-act_ori)/eps;
    double nu_minus=-(act_ori-act_minus)/eps;
    double nu=-(act_plus-act_minus)/(2*eps);
    
    vector_reset(F);
#endif
    
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
		
		su3_summ_the_prod_double(F[ivol][mu],E,-kappa/8); //8!?
	      }
	  }
	THREAD_BARRIER();
      }
    set_borders_invalid(F);
    
#ifdef DEBUG
    su3 r1,r2;
    unsafe_su3_prod_su3(r1,conf[0][0],F[0][0]);
    unsafe_su3_traceless_anti_hermitian_part(r2,r1);
    double tr=(r2[1][0][1]+r2[0][1][1])/2;
    printf("an: %+016.016le, nu: %+016.016le, nu+: %+016.016le, nu-: %+016.016le\n",tr,nu,nu_plus,nu_minus);
    crash("anna QCD force originating from H");
#endif
    
    nissa_free(temp);
    nissa_free(H_nu);
  }
  THREADABLE_FUNCTION_END
}
