#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/edges.hpp"
#include "geometry/geometry_eo.hpp"
#include "linalgs/reduce.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //compute plaquettes and rectangles
  THREADABLE_FUNCTION_2ARG(global_plaquette_and_rectangles_eo_conf, double*,glb_shapes, eo_ptr<quad_su3>,conf)
  {
    GET_THREAD_ID();
    
    communicate_eo_quad_su3_edges(conf);
    
    //summ squares and rectangles separately
    complex *point_shapes=nissa_malloc("point_shapes",loc_vol,complex);
    vector_reset(point_shapes);
    
    for(int par=0;par<2;par++)
      for(int mu=0;mu<NDIM;mu++) //link dir
	for(int nu=0;nu<NDIM;nu++) //staple dir
	  if(nu!=mu)
	    {
	      NISSA_PARALLEL_LOOP(A,0,loc_volh)
		{
		  int ivol=loclx_of_loceo[par][A];
		  
		  //compute forward staple starting from A
		  int B=loceo_neighup[par][A][nu],D=loceo_neighdw[par][A][nu];
		  int E=loceo_neighup[!par][D][mu],F=loceo_neighup[par][A][mu];
		  su3 ABC,ABCF;
		  unsafe_su3_prod_su3(ABC,conf[par][A][nu],conf[!par][B][mu]);
		  unsafe_su3_prod_su3_dag(ABCF,ABC,conf[!par][F][nu]);
		  
		  //taking the trace we summ to plaq_summ (only if nu>mu)
		  if(nu>mu) point_shapes[ivol][RE]+=real_part_of_trace_su3_prod_su3_dag(ABCF,conf[par][A][mu]);
		  
		  //compute backward staple starting from A
		  su3 ADE,ADEF;
		  unsafe_su3_dag_prod_su3(ADE,conf[!par][D][nu],conf[!par][D][mu]);
		  unsafe_su3_prod_su3(ADEF,ADE,conf[par][E][nu]);
		  
		  //taking the trace we summ to rect_summ
		  point_shapes[ivol][IM]+=real_part_of_trace_su3_prod_su3_dag(ABCF,ADEF);
		}
	      NISSA_PARALLEL_LOOP_END;
	    }
    THREAD_BARRIER();
    
    //reduce and free
    complex coll_shapes;
    glb_reduce(&coll_shapes,point_shapes,loc_vol);
    nissa_free(point_shapes);
    
    //normalize (passing throug additional var because of external unkwnon env)
    glb_shapes[RE]=coll_shapes[RE]/(18*glb_vol);
    glb_shapes[IM]=coll_shapes[IM]/(36*glb_vol);
  }
  THREADABLE_FUNCTION_END

  //compute plaquettes and rectangles
  THREADABLE_FUNCTION_2ARG(point_plaquette_and_rectangles_lx_conf, complex*,point_shapes, quad_su3*,conf)
  {
    GET_THREAD_ID();
    
    //communicate conf and reset point shapes
    communicate_lx_quad_su3_edges(conf);
    vector_reset(point_shapes);
    
    for(int mu=0;mu<4;mu++) //link dir
      for(int nu=0;nu<4;nu++) //staple dir
	if(nu!=mu)
	  {
	    NISSA_PARALLEL_LOOP(A,0,loc_vol)
	      {
		int ivol=A;
		
		//compute forward staple starting from A
		int B=loclx_neighup[A][nu],D=loclx_neighdw[A][nu];
		int E=loclx_neighup[D][mu],F=loclx_neighup[A][mu];
		su3 ABC,ABCF;
		unsafe_su3_prod_su3(ABC,conf[A][nu],conf[B][mu]);
		unsafe_su3_prod_su3_dag(ABCF,ABC,conf[F][nu]);
		
		//taking the trace we summ to plaq_summ (only if nu>mu)
		if(nu>mu) point_shapes[ivol][RE]+=real_part_of_trace_su3_prod_su3_dag(ABCF,conf[A][mu]);
		
		//compute backward staple starting from A
		su3 ADE,ADEF;
		unsafe_su3_dag_prod_su3(ADE,conf[D][nu],conf[D][mu]);
		unsafe_su3_prod_su3(ADEF,ADE,conf[E][nu]);
		
		//taking the trace we summ to rect_summ
		point_shapes[ivol][IM]+=real_part_of_trace_su3_prod_su3_dag(ABCF,ADEF);
	      }
	    NISSA_PARALLEL_LOOP_END;
	  }
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
  
  //compute plaquettes and rectangles
  THREADABLE_FUNCTION_2ARG(global_plaquette_and_rectangles_lx_conf, double*,glb_shapes, quad_su3*,conf)
  {
    //summ squares and rectangles separately
    complex *point_shapes=nissa_malloc("point_shapes",loc_vol,complex);
    point_plaquette_and_rectangles_lx_conf(point_shapes,conf);
    
    //reduce and free
    complex coll_shapes;
    glb_reduce(&coll_shapes,point_shapes,loc_vol);
    nissa_free(point_shapes);
    
    //normalize (passing throug additional var because of external unkwnon env)
    glb_shapes[RE]=coll_shapes[RE]/(18*glb_vol);
    glb_shapes[IM]=coll_shapes[IM]/(36*glb_vol);
  }
  THREADABLE_FUNCTION_END
  
  //compute plaquettes and rectangles
  THREADABLE_FUNCTION_2ARG(global_plaquette_and_rectangles_lx_conf_per_timeslice, double*,glb_shapes, quad_su3*,conf)
  {
    GET_THREAD_ID();
    
    //summ squares and rectangles separately
    complex *point_shapes=nissa_malloc("point_shapes",loc_vol,complex);
    point_plaquette_and_rectangles_lx_conf(point_shapes,conf);
    
    //reduce
    complex *loc_shapes=nissa_malloc("loc_shapes",glb_size[0],complex);
    vector_reset(loc_shapes);
    
    //loop over time
    NISSA_PARALLEL_LOOP(loc_t,0,loc_size[0])
      for(int ivol=loc_t*loc_spat_vol;ivol<(loc_t+1)*loc_spat_vol;ivol++)
	complex_summassign(loc_shapes[glb_coord_of_loclx[ivol][0]],point_shapes[ivol]);
    NISSA_PARALLEL_LOOP_END;
    nissa_free(point_shapes);
    
    //reduce (passing throug additional var because of external unkwnon env)
    complex *coll_shapes=nissa_malloc("coll_shapes",glb_size[0],complex);
    if(IS_MASTER_THREAD) MPI_Reduce(loc_shapes,coll_shapes,2*glb_size[0],MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    nissa_free(loc_shapes);
    
    //normalize
    for(int t=0;t<glb_size[0];t++)
      {
	glb_shapes[2*t+0]=coll_shapes[t][RE]/(18*glb_vol/glb_size[0]);
	glb_shapes[2*t+1]=coll_shapes[t][IM]/(36*glb_vol/glb_size[0]);
      }
    nissa_free(coll_shapes);
  }
  THREADABLE_FUNCTION_END
}
