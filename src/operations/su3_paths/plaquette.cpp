#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/dirac.hpp"
#include "new_types/float_128.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3_op.hpp"
#include "routines/mpi_routines.hpp"

namespace nissa
{
  /*
    
      C------D
    n |      |
    u |      |
      A--mu--B
      
      The square path P_{mu,nu} is defined as U(A,mu)U(B,nu)U^(C,mu)U^(A,nu)=
      =U(A,mu)U(B,nu)(U(A,nu)U^(C,mu))^=U(AB,munu)*U^(AC,numu)
  */
  
  CUDA_HOST_AND_DEVICE void point_plaquette_lx_conf(complex loc_plaq,quad_su3 *conf,int A)
  {
    loc_plaq[0]=loc_plaq[1]=0;
    for(int mu=0;mu<NDIM;mu++)
      {
	int B=loclx_neighup[A][mu];
	for(int nu=mu+1;nu<NDIM;nu++)
	  {
	    int C=loclx_neighup[A][nu];
	    su3 ABD,ACD;
	    unsafe_su3_prod_su3(ABD,conf[A][mu],conf[B][nu]);
	    unsafe_su3_prod_su3(ACD,conf[A][nu],conf[C][mu]);
	    
	    int ts=(mu!=0&&nu!=0);
	    loc_plaq[ts]+=real_part_of_trace_su3_prod_su3_dag(ABD,ACD);
	  }
      }
  }
  CUDA_HOST_AND_DEVICE void point_plaquette_eo_conf(complex loc_plaq,eo_ptr<quad_su3> conf,int par,int A)
  {
    loc_plaq[0]=loc_plaq[1]=0;
    for(int mu=0;mu<NDIM;mu++)
      {
	int B=loceo_neighup[par][A][mu];
	for(int nu=mu+1;nu<NDIM;nu++)
	  {
	    int C=loceo_neighup[par][A][nu];
	    su3 ABD,ACD;
	    unsafe_su3_prod_su3(ABD,conf[par][A][mu],conf[!par][B][nu]);
	    unsafe_su3_prod_su3(ACD,conf[par][A][nu],conf[!par][C][mu]);
	    
	    int ts=(mu!=0&&nu!=0);
	    loc_plaq[ts]+=real_part_of_trace_su3_prod_su3_dag(ABD,ACD);
	  }
      }
  }
  
  //calculate the global plaquette of an lx conf
  THREADABLE_FUNCTION_2ARG(global_plaquette_lx_conf, double*,totplaq, quad_su3*,conf)
  {
    GET_THREAD_ID();
    
    //summ temporal and spatial separately
    complex *point_plaq=nissa_malloc("point_plaq",loc_vol,complex);
    communicate_lx_quad_su3_borders(conf);
    
    //loop over all the lattice
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      point_plaquette_lx_conf(point_plaq[ivol],conf,ivol);
    NISSA_PARALLEL_LOOP_END;
    
    //wait to have filled all the point array
    THREAD_BARRIER();
    
    //reduce as complex and normalize
    complex temp;
    complex_vector_glb_collapse(temp,point_plaq,loc_vol);
    for(int ts=0;ts<2;ts++) totplaq[ts]=temp[ts]/(glb_vol*3*3);
    
    nissa_free(point_plaq);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_2ARG(global_plaquette_eo_conf, double*,totplaq, eo_ptr<quad_su3>,conf)
  {
    GET_THREAD_ID();
    
    //summ temporal and spatial separately
    complex *point_plaq=nissa_malloc("point_plaq",loc_vol,complex);
    #warning communicate_ev_and_od_quad_su3_borders(conf);
    
    //loop over all the lattice
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  point_plaquette_eo_conf(point_plaq[loclx_of_loceo[par][ieo]],conf,par,ieo);
	NISSA_PARALLEL_LOOP_END;
      }
    
    //wait to have filled all the point array
    THREAD_BARRIER();
    
    //reduce as complex and normalize
    complex temp;
    complex_vector_glb_collapse(temp,point_plaq,loc_vol);
    for(int ts=0;ts<2;ts++) totplaq[ts]=temp[ts]/(glb_vol*3*3);
    
    nissa_free(point_plaq);
  }
  THREADABLE_FUNCTION_END
  
  //return the average between spatial and temporary plaquette
  double global_plaquette_lx_conf(quad_su3 *conf)
  {
    double plaq[2];
    global_plaquette_lx_conf(plaq,conf);
    return (plaq[0]+plaq[1])/2;
  }
  double global_plaquette_eo_conf(eo_ptr<quad_su3> conf)
  {
    //compute the two plaquettes
    double plaq[2];
    global_plaquette_eo_conf(plaq,conf);
    return (plaq[0]+plaq[1])/2;
  }
}
