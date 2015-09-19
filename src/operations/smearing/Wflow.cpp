#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "communicate/communicate.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //flow for the given time for a dt using 1006.4518 appendix C
  THREADABLE_FUNCTION_2ARG(Wflow_lx_conf, quad_su3*,conf, double,dt)
  {
    GET_THREAD_ID();
    
    //storage for staples
    quad_su3 *arg=nissa_malloc("arg",loc_vol,quad_su3);
    vector_reset(arg);
    
    //we write the 4 terms of the Runge Kutta scheme iteratively
    //we add with the new weight the previous one multiplied by the old weight
    double RK_wn[3]={1.0/4,8.0/9, 3.0/4};
    double RK_wo[3]={0,    -17/9, -1};
    
    for(int iter=0;iter<3;iter++)
      {
	communicate_lx_quad_su3_edges(conf);
	
	//add the new argument of the exponential to the old one
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      //compute the new contribution
	      su3 staple,temp;
	      su3_put_to_zero(staple);
	      for(int inu=0;inu<NDIM-1;inu++)
		{
		  int nu=perp_dir[mu][inu];
		  int A=ivol,B=loclx_neighup[A][nu],D=loclx_neighdw[A][nu],E=loclx_neighup[D][mu],F=loclx_neighup[A][mu];
		  unsafe_su3_prod_su3(       temp, conf[A][nu],conf[B][mu]);
		  su3_summ_the_prod_su3_dag(staple,temp,       conf[F][nu]);
		  unsafe_su3_dag_prod_su3(temp,    conf[D][nu],conf[D][mu]);
		  su3_summ_the_prod_su3(staple,    temp,       conf[E][nu]);
		}
	      
	      //build Omega
	      su3 omega;
	      unsafe_su3_prod_su3_dag(omega,staple,conf[ivol][mu]);
	      
	      //compute Q and weight (the minus is there due to original stout
	      su3 iQ,Q;
	      unsafe_su3_traceless_anti_hermitian_part(iQ,omega);
	      su3_prod_idouble(Q,iQ,-RK_wn[iter]*dt); //putting here the integration time
	      
	      //combine old and new
	      su3_prod_double(arg[ivol][mu],arg[ivol][mu],RK_wo[iter]);
	      su3_summassign(arg[ivol][mu],Q);
	    }
	THREAD_BARRIER();
	
	//integrate
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      su3 expiQ;
	      safe_anti_hermitian_exact_i_exponentiate(expiQ,arg[ivol][mu]);
	      safe_su3_prod_su3(conf[ivol][mu],expiQ,conf[ivol][mu]);
	    }
	set_borders_invalid(conf);
      }
    
    nissa_free(arg);
  }
  THREADABLE_FUNCTION_END
  
  //flow for the given time using an adaptative integrator of first order
  THREADABLE_FUNCTION_4ARG(adaptative_stout_lx_conf, quad_su3*,conf, double*,t, double,Tmax, double,arg_max)
  {
    GET_THREAD_ID();
    
    //storage for staples
    quad_su3 *arg=nissa_malloc("arg",loc_vol,quad_su3);
    vector_reset(arg);
    
    communicate_lx_quad_su3_edges(conf);
    
    //add the new argument of the exponential to the old one
    double force_norm=0;
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  //compute the new contribution
	  su3 staple,temp;
	  su3_put_to_zero(staple);
	  for(int inu=0;inu<NDIM-1;inu++)
	    {
	      int nu=perp_dir[mu][inu];
	      int A=ivol,B=loclx_neighup[A][nu],D=loclx_neighdw[A][nu],E=loclx_neighup[D][mu],F=loclx_neighup[A][mu];
	      unsafe_su3_prod_su3(       temp, conf[A][nu],conf[B][mu]);
	      su3_summ_the_prod_su3_dag(staple,temp,       conf[F][nu]);
	      unsafe_su3_dag_prod_su3(temp,    conf[D][nu],conf[D][mu]);
	      su3_summ_the_prod_su3(staple,    temp,       conf[E][nu]);
	    }
	  
	  //build Omega
	  su3 omega;
	  unsafe_su3_prod_su3_dag(omega,staple,conf[ivol][mu]);
	  
	  //compute Q and weight (the minus is there due to original stout
	  unsafe_su3_traceless_anti_hermitian_part(arg[ivol][mu],omega);
	  force_norm+=su3_norm2(arg[ivol][mu]);
	}
    force_norm=sqrt(glb_reduce_double(force_norm)/(NDIM*glb_vol));
    double dt=std::min((Tmax-*t),arg_max/force_norm);
    master_printf("Force norm: %lg, dt: %lg\n",force_norm,dt);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  su3 Q;
	  su3_prod_idouble(Q,arg[ivol][mu],-dt); //putting here the integration time (- due to i)
	}
    THREAD_BARRIER();
    
    //integrate
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  su3 expiQ;
	  safe_anti_hermitian_exact_i_exponentiate(expiQ,arg[ivol][mu]);
	  safe_su3_prod_su3(conf[ivol][mu],expiQ,conf[ivol][mu]);
	}
    
    if(IS_MASTER_THREAD) *t+=dt;
    set_borders_invalid(conf);
    
    nissa_free(arg);
  }
  THREADABLE_FUNCTION_END

}
