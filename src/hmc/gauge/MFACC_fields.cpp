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
    double glb_action_lx[2];
    for(int ifield=0;ifield<2;ifield++)
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
    
    for(int ifield=0;ifield<2;ifield++)
      {
        //compute
        apply_MFACC(temp,conf,kappa,pi[ifield]);
        apply_MFACC(F,conf,kappa,temp);
        
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
  
  THREADABLE_FUNCTION_4ARG(compa, su3*,out, quad_su3*,conf, double,kappa, su3*,in)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
        //reset
        su3_put_to_zero(out[ivol]);
        
        for(int mu=0;mu<NDIM;mu++)
          {
            //neighbours search
            int iup=loclx_neighup[ivol][mu];
            int idw=loclx_neighdw[ivol][mu];
            
            su3 temp;
            
	    //su3_summ_the_prod_su3(out[ivol],conf[ivol][mu],in[iup]);
            unsafe_su3_prod_su3(temp,conf[ivol][mu],in[iup]);
            su3_summ_the_prod_su3_dag(out[ivol],temp,conf[ivol][mu]);
            
            //su3_summ_the_dag_prod_su3(out[ivol],conf[idw][mu],in[idw]);
            unsafe_su3_dag_prod_su3(temp,conf[idw][mu],in[idw]);
            su3_summ_the_prod_su3(out[ivol],temp,conf[idw][mu]);
          }
      }
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //compute the MFACC momenta-related QCD force (derivative of \pi^\dag MM \pi / 2 )
  THREADABLE_FUNCTION_5ARG(MFACC_momenta_QCD_force, quad_su3*,F, quad_su3*,conf, double,kappa, su3**,pi, bool,reset)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Computing Fourier QCD force originated by acceleration-momenta\n");
    
    su3 *temp=nissa_malloc("temp",loc_vol+bord_vol,su3);
    if(reset) vector_reset(F);
    
#ifdef DEBUG
    vector_reset(F);
    
    double eps=1e-5;
    
    //random gauge transform
    su3 *fixm=nissa_malloc("fixm",loc_vol+bord_vol,su3);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      su3_put_to_id(fixm[ivol]);
      //su3_put_to_rnd(fixm[ivol],loc_rnd_gen[ivol]);
    set_borders_invalid(fixm);
    gauge_transform_conf(conf,fixm,conf);
    for(int ifield=0;ifield<2;ifield++)
      NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	{
	  //safe_su3_prod_su3(pi[ifield][ivol],fixm[ivol],pi[ifield][ivol]);
	  su3 a; //gauge transformation when more true operator used
	  unsafe_su3_prod_su3(a,fixm[ivol],pi[ifield][ivol]);
	  unsafe_su3_prod_su3_dag(pi[ifield][ivol],a,fixm[ivol]);
	}
    
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
    
    for(int ifield=0;ifield<2;ifield++)
      {
	apply_MFACC(temp,conf,kappa,pi[ifield]);
	//compa(temp,conf,kappa,pi[ifield]);
	
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      su3 t;
	      int up=loclx_neighup[ivol][mu];
	      
	      su3 E;
	      
	      //forward piece
	      unsafe_su3_dag_prod_su3_dag(t,conf[ivol][mu],temp[ivol]);
	      unsafe_su3_prod_su3(E,pi[ifield][up],t);
	      
	      //backward piece
	      unsafe_su3_dag_prod_su3_dag(t,temp[up],conf[ivol][mu]);
	      su3_summ_the_prod_su3(E,t,pi[ifield][ivol]);
	      
	      su3_summ_the_prod_double(F[ivol][mu],E,kappa/8); //factor of 2
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
    nissa_free(fixm);
crash("anna");
#endif
  }
  THREADABLE_FUNCTION_END
  
  //compute the MFACC fields-related QCD force (derivative of \H^\dag G \H/2)
  THREADABLE_FUNCTION_7ARG(MFACC_QCD_momenta_QCD_force, quad_su3*,F, quad_su3*,conf, double,kappa, int,niter, double,residue, quad_su3*,H, bool,reset)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Computing Fourier acceleration fields originated QCD force\n");
    
    su3 *H_nu=nissa_malloc("H_nu",loc_vol+bord_vol,su3);
    su3 *temp=nissa_malloc("temp",loc_vol+bord_vol,su3);
    if(reset) vector_reset(F);
    
#ifdef DEBUG
    double eps=1.e-5;
    
    //store original
    su3 sto;
    su3_copy(sto,conf[0][0]);
    
    //compute action previous to change
    double pre=0;
    inv_MFACC_cg(temp,NULL,conf,kappa,niter,residue,H_nu);
    double_vector_glb_scalar_prod(&(pre),(double*)H_nu,(double*)temp,18*loc_vol);
    
    //prepare increment and change
    su3 mod,ba;
    su3_put_to_zero(ba);
    ba[1][0][0]=ba[0][1][0]=eps/4;
    safe_anti_hermitian_exact_i_exponentiate(mod,ba);
    safe_su3_prod_su3(conf[0][0],mod,sto);
    
    //compute it after
    double post=0;
    inv_MFACC_cg(temp,NULL,conf,kappa,niter,residue,H_nu);
    double_vector_glb_scalar_prod(&(post),(double*)H_nu,(double*)temp,18*loc_vol);
    
    //set back everything
    printf("pre: %lg, post: %lg\n",pre,post);
    su3_copy(conf[0][0],sto);
    //double num=-(post-pre)/eps;
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
	    su3 E;
	    su3_put_to_zero(E);
	    
	    for(int mu=0;mu<NDIM;mu++)
	      {
		su3 t;
		int up=loclx_neighup[ivol][mu];
		
		//version contracting with temp
		//forward piece
		unsafe_su3_dag_prod_su3_dag(t,conf[ivol][mu],temp[ivol]);
		su3_summ_the_prod_su3(E,temp[up],t);
		unsafe_su3_dag_prod_su3(t,conf[ivol][mu],temp[ivol]);
		su3_summ_the_dag_prod_su3(E,temp[up],t);
		//backward piece
		unsafe_su3_dag_prod_su3_dag(t,temp[up],conf[ivol][mu]);
		su3_summ_the_prod_su3(E,t,temp[ivol]);
		unsafe_su3_prod_su3_dag(t,temp[up],conf[ivol][mu]);
		su3_summ_the_prod_su3_dag(E,t,temp[ivol]);
	      }
	    
	    su3_summ_the_prod_double(F[ivol][nu],E,-kappa/32);
	    
	    //version contracting with pi
	    //forward piece
	    //unsafe_su3_dag_prod_su3_dag(t,conf[ivol][mu],pi[ifield][ivol]);
	    //su3_summ_the_prod_su3(F[ivol][mu],pi[ifield][up],t);
	    //backward piece
	    //unsafe_su3_prod_su3_dag(t,pi[ifield][up],conf[ivol][mu]);
	    //su3_summ_the_prod_su3_dag(F[ivol][mu],t,pi[ifield][ivol]);
	    
	    //su3_summ_the_prod_su3_dag(F[ivol][mu],pi[0][loclx_neighup[ivol][mu]],pi[0][ivol]);
	  }
	THREAD_BARRIER();
      }
    set_borders_invalid(F);
    
    nissa_free(temp);
    nissa_free(H_nu);
    
    /*
    su3 r1,r2;
    unsafe_su3_prod_su3(r1,conf[0][0],F[0][0]);
    unsafe_su3_traceless_anti_hermitian_part(r2,r1);
    double tr=(r2[1][0][1]+r2[0][1][1])/4;
    printf("an: %lg, nu: %lg\n",tr,nu);
    nissa_free(fixm);
    */
  }
  THREADABLE_FUNCTION_END
}
