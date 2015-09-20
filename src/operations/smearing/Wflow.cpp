#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "communicate/communicate.hpp"
#include "new_types/su3.hpp"
#include "operations/su3_paths/plaquette.hpp"
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
  THREADABLE_FUNCTION_4ARG(adaptative_stout_lx_conf, quad_su3*,conf, double*,t, double,Tmax, double*,ext_dt)
  {
    GET_THREAD_ID();
    
    //storage for staples
    quad_su3 *arg=nissa_malloc("arg",loc_vol,quad_su3);
    quad_su3 *test_conf=nissa_malloc("test_conf",loc_vol,quad_su3);
    vector_reset(arg);
    
    communicate_lx_quad_su3_edges(conf);
    
    //add the new argument of the exponential to the old one
    double plaq=0;
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
	  plaq+=su3_real_trace(omega)/(2*NDIM-2);
	  
	  //compute Q and weight (the minus is there due to original stout
	  unsafe_su3_traceless_anti_hermitian_part(arg[ivol][mu],omega);
	}
    plaq=glb_reduce_double(plaq)/NDIM/NCOL/glb_vol;
    
    //dt=0.10;
    
    double dt=0;
    double tstep=*ext_dt;
    double ttoll=0.001;
    double old_plaq=global_plaquette_lx_conf(conf);
    master_printf("Plaq: %16.16lg, at t=0\n",old_plaq);
    
    do
      {
	double dt_test=dt+tstep;
	
	//integrate
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      su3 Q,expiQ;
	      su3_prod_idouble(Q,arg[ivol][mu],-dt_test);
	      safe_anti_hermitian_exact_i_exponentiate(expiQ,Q);
	      safe_su3_prod_su3(test_conf[ivol][mu],expiQ,conf[ivol][mu]);
	    }
	set_borders_invalid(test_conf);
	double new_plaq=global_plaquette_lx_conf(test_conf);
        master_printf("Plaq: %16.16lg vs %16.16lg, at t=%d\n",new_plaq,old_plaq,dt_test);
	
	if(new_plaq>old_plaq)
	  {
	    old_plaq=new_plaq;
	    dt=dt_test;
	    master_printf(" accepted\n");
	  }
	else
	  {
	    tstep/=2;
	    master_printf(" rejected\n");
	  }
      }
    while(tstep>=ttoll);
    
    nissa_free(test_conf);
    
    //integrate
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  su3 Q,expiQ;
	  su3_prod_idouble(Q,arg[ivol][mu],-dt);
	  safe_anti_hermitian_exact_i_exponentiate(expiQ,Q);
	  safe_su3_prod_su3(conf[ivol][mu],expiQ,conf[ivol][mu]);
	}
    
    if(IS_MASTER_THREAD)
      {
	*t+=dt;
	*ext_dt=dt;
      }
    set_borders_invalid(conf);
    
    nissa_free(arg);
  }
  THREADABLE_FUNCTION_END
}
