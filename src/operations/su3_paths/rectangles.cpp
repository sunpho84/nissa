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
  void global_plaquette_and_rectangles_eo_conf(double* glb_shapes,eo_ptr<quad_su3> conf)
  {
    
    communicate_eo_quad_su3_edges(conf);
    
    //summ squares and rectangles separately
    complex *point_shapes=nissa_malloc("point_shapes",locVol.nastyConvert(),complex);
    vector_reset(point_shapes);
    
    for(int par=0;par<2;par++)
      for(int mu=0;mu<NDIM;mu++) //link dir
	for(int nu=0;nu<NDIM;nu++) //staple dir
	  if(nu!=mu)
	    {
	      NISSA_PARALLEL_LOOP(A,0,locVolh)
		{
		  const LocLxSite ivol=loclx_of_loceo[par][A];
		  
		  //compute forward staple starting from A
		  int B=loceo_neighup[par][A][nu],D=loceo_neighdw[par][A][nu];
		  int E=loceo_neighup[!par][D][mu],F=loceo_neighup[par][A][mu];
		  su3 ABC,ABCF;
		  unsafe_su3_prod_su3(ABC,conf[par][A][nu],conf[!par][B][mu]);
		  unsafe_su3_prod_su3_dag(ABCF,ABC,conf[!par][F][nu]);
		  
		  //taking the trace we summ to plaq_summ (only if nu>mu)
		  if(nu>mu) point_shapes[ivol.nastyConvert()][RE]+=real_part_of_trace_su3_prod_su3_dag(ABCF,conf[par][A][mu]);
		  
		  //compute backward staple starting from A
		  su3 ADE,ADEF;
		  unsafe_su3_dag_prod_su3(ADE,conf[!par][D][nu],conf[!par][D][mu]);
		  unsafe_su3_prod_su3(ADEF,ADE,conf[par][E][nu]);
		  
		  //taking the trace we summ to rect_summ
		  point_shapes[ivol.nastyConvert()][IM]+=real_part_of_trace_su3_prod_su3_dag(ABCF,ADEF);
		}
	      NISSA_PARALLEL_LOOP_END;
	    }
    THREAD_BARRIER();
    
    //reduce and free
    complex coll_shapes;
    glb_reduce(&coll_shapes,point_shapes,locVol.nastyConvert());
    nissa_free(point_shapes);
    
    //normalize (passing throug additional var because of external unkwnon env)
    glb_shapes[RE]=coll_shapes[RE]/(18*glbVol());
    glb_shapes[IM]=coll_shapes[IM]/(36*glbVol());
  }

  //compute plaquettes and rectangles
  void point_plaquette_and_rectangles_lx_conf(complex* point_shapes,quad_su3* conf)
  {
    
    //communicate conf and reset point shapes
    communicate_lx_quad_su3_edges(conf);
    vector_reset(point_shapes);
    
    for(int mu=0;mu<4;mu++) //link dir
      for(int nu=0;nu<4;nu++) //staple dir
	if(nu!=mu)
	  {
	    NISSA_PARALLEL_LOOP(A,0,locVol)
	      {
		const LocLxSite ivol=A;
		
		//compute forward staple starting from A
		int B=loclxNeighup[A.nastyConvert()][nu],D=loclxNeighdw[A.nastyConvert()][nu];
		int E=loclxNeighup[D][mu],F=loclxNeighup[A.nastyConvert()][mu];
		su3 ABC,ABCF;
		unsafe_su3_prod_su3(ABC,conf[A.nastyConvert()][nu],conf[B][mu]);
		unsafe_su3_prod_su3_dag(ABCF,ABC,conf[F][nu]);
		
		//taking the trace we summ to plaq_summ (only if nu>mu)
		if(nu>mu) point_shapes[ivol.nastyConvert()][RE]+=real_part_of_trace_su3_prod_su3_dag(ABCF,conf[A.nastyConvert()][mu]);
		
		//compute backward staple starting from A
		su3 ADE,ADEF;
		unsafe_su3_dag_prod_su3(ADE,conf[D][nu],conf[D][mu]);
		unsafe_su3_prod_su3(ADEF,ADE,conf[E][nu]);
		
		//taking the trace we summ to rect_summ
		point_shapes[ivol.nastyConvert()][IM]+=real_part_of_trace_su3_prod_su3_dag(ABCF,ADEF);
	      }
	    NISSA_PARALLEL_LOOP_END;
	  }
    THREAD_BARRIER();
  }
  
  //compute plaquettes and rectangles
  void global_plaquette_and_rectangles_lx_conf(double* glb_shapes,quad_su3* conf)
  {
    //summ squares and rectangles separately
    complex *point_shapes=nissa_malloc("point_shapes",locVol.nastyConvert(),complex);
    point_plaquette_and_rectangles_lx_conf(point_shapes,conf);
    
    //reduce and free
    complex coll_shapes;
    glb_reduce(&coll_shapes,point_shapes,locVol.nastyConvert());
    nissa_free(point_shapes);
    
    //normalize (passing throug additional var because of external unkwnon env)
    glb_shapes[RE]=coll_shapes[RE]/(18*glbVol());
    glb_shapes[IM]=coll_shapes[IM]/(36*glbVol());
  }
  
  //compute plaquettes and rectangles
  void global_plaquette_and_rectangles_lx_conf_per_timeslice(double* glb_shapes,quad_su3* conf)
  {
    
    //summ squares and rectangles separately
    complex *point_shapes=nissa_malloc("point_shapes",locVol.nastyConvert(),complex);
    point_plaquette_and_rectangles_lx_conf(point_shapes,conf);
    
    //reduce
    complex *loc_shapes=nissa_malloc("loc_shapes",glbSize[0],complex);
    vector_reset(loc_shapes);
    
    //loop over time
    NISSA_PARALLEL_LOOP(loc_t,0,locSize[0])
      for(LocLxSite ivol=loc_t*locSpatVol;ivol<(loc_t+1)*locSpatVol;ivol++)
	complex_summassign(loc_shapes[glbCoordOfLoclx[ivol.nastyConvert()][0]],point_shapes[ivol.nastyConvert()]);
    NISSA_PARALLEL_LOOP_END;
    nissa_free(point_shapes);
    
    //reduce (passing throug additional var because of external unkwnon env)
    complex *coll_shapes=nissa_malloc("coll_shapes",glbSize[0],complex);
    if(IS_MASTER_THREAD) MPI_Reduce(loc_shapes,coll_shapes,2*glbSize[0],MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    nissa_free(loc_shapes);
    
    //normalize
    for(int t=0;t<glbSize[0];t++)
      {
	glb_shapes[2*t+0]=coll_shapes[t][RE]/(18.0*glbVol()/glbSize[0]);
	glb_shapes[2*t+1]=coll_shapes[t][IM]/(36.0*glbVol()/glbSize[0]);
      }
    nissa_free(coll_shapes);
  }
}
