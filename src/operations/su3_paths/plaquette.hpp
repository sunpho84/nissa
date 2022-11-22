#ifndef _PLAQUETTE_HPP
#define _PLAQUETTE_HPP

#include "base/field.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "linalgs/reduce.hpp"
#include "new_types/su3_op.hpp"

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
  template <typename Plaq,
	    typename Conf>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void point_plaquette_lx_conf(Plaq&& loc_plaq,
			       const Conf& conf,
			       const int& A)
  {
    loc_plaq[0]=loc_plaq[1]=0.0;
    
    for(int mu=0;mu<NDIM;mu++)
      {
	const int B=loclxNeighup[A][mu];
	
	for(int nu=mu+1;nu<NDIM;nu++)
	  {
	    const int C=loclxNeighup[A][nu];
	    
	    su3 ABD;
	    unsafe_su3_prod_su3(ABD,conf[A][mu],conf[B][nu]);
	    
	    su3 ACD;
	    unsafe_su3_prod_su3(ACD,conf[A][nu],conf[C][mu]);
	    
	    const int ts=(mu!=0 and nu!=0);
	    loc_plaq[ts]+=real_part_of_trace_su3_prod_su3_dag(ABD,ACD);
	  }
      }
  }
  
  /// Computes the global plaquette, separating timelike and not
  template <typename C>
  void global_plaquette_lx_conf(double* totplaq,
				const C& conf)
  {
    if constexpr(std::is_pointer_v<C>) // hack
      communicate_lx_quad_su3_borders(conf);
    else
      conf.updateHalo();
    
    Field<complex> point_plaq("point_plaq");
    
    //loop over all the lattice
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      point_plaquette_lx_conf(point_plaq[ivol],conf,ivol);
    NISSA_PARALLEL_LOOP_END;
    
    //reduce as complex and normalize
    complex temp;
    glb_reduce(&temp,point_plaq,locVol);
    for(int ts=0;ts<2;ts++)
      totplaq[ts]=temp[ts]/(glbVol*NCOL*NCOL);
  }
  
  /// Average plaquette
  template <typename C>
  double global_plaquette_lx_conf(const C& conf)
    {
      double plaq[2];
      
      global_plaquette_lx_conf(plaq,conf);
      
      return (plaq[0]+plaq[1])/2;
    }
  
  double global_plaquette_eo_conf(eo_ptr<quad_su3> conf);
  void global_plaquette_eo_conf(double *totplaq,eo_ptr<quad_su3> conf);
}

#endif
