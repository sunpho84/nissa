#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //apply DD
  THREADABLE_FUNCTION_5ARG(apply_MFACC, quad_su3*,out, quad_su3*,conf, double,kappa, double,offset, quad_su3*,in)
  {
    if(!check_borders_valid(in)) communicate_lx_quad_su3_borders(in);
    if(!check_borders_valid(conf)) communicate_lx_quad_su3_borders(conf);
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	for(int nu=0;nu<NDIM;nu++) su3_put_to_zero(out[ivol][nu]);
	
	for(int mu=0;mu<NDIM;mu++)
	  {
	    //neighbours search
	    int iup=loclx_neighup[ivol][mu];
	    int idw=loclx_neighdw[ivol][mu];
	    
	    for(int nu=0;nu<NDIM;nu++)
	      {
		su3 temp;
		
		unsafe_su3_prod_su3(temp,conf[ivol][mu],in[iup][nu]);
		su3_summ_the_prod_su3_dag(out[ivol][nu],temp,conf[ivol][mu]);
		
		unsafe_su3_dag_prod_su3(temp,conf[idw][mu],in[idw][nu]);
		su3_summ_the_prod_su3(out[ivol][nu],temp,conf[idw][mu]);
	      }
	  }
	
	//add normalization and similar
	for(int nu=0;nu<NDIM;nu++)
	  {
	    su3_prod_double(out[ivol][nu],out[ivol][nu],-kappa/(4*NDIM));
	    su3_summ_the_prod_double(out[ivol][nu],in[ivol][nu],1-kappa/2+offset);
	  }
      }
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_5ARG(apply_MFACC, su3*,out, quad_su3*,conf, double,kappa, double,offset, su3*,in)
  {
    if(!check_borders_valid(in)) communicate_lx_su3_borders(in);
    if(!check_borders_valid(conf)) communicate_lx_quad_su3_borders(conf);
    
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
	    
	    unsafe_su3_prod_su3(temp,conf[ivol][mu],in[iup]);
	    su3_summ_the_prod_su3_dag(out[ivol],temp,conf[ivol][mu]);
	    
	    unsafe_su3_dag_prod_su3(temp,conf[idw][mu],in[idw]);
	    su3_summ_the_prod_su3(out[ivol],temp,conf[idw][mu]);
	  }
	
	//add normalization and similar
	su3_prod_double(out[ivol],out[ivol],-kappa/(4*NDIM));
	su3_summ_the_prod_double(out[ivol],in[ivol],1-kappa/2+offset);
      }
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
}
