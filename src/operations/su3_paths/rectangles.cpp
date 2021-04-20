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
    
    FOR_BOTH_PARITIES(par)
      FOR_ALL_DIRS(mu) //link dir
	FOR_ALL_DIRS(nu) //staple dir
	  if(nu!=mu)
	    {
	      NISSA_PARALLEL_LOOP(A,0,locVolh)
		{
		  const LocLxSite& ivol=loclx_of_loceo(par,A);
		  
		  //compute forward staple starting from A
		  const LocEoSite& B=loceo_neighup(par,A,nu),D=loceo_neighdw(par,A,nu);
		  const LocEoSite& E=loceo_neighup(1-par,D,mu),F=loceo_neighup(par,A,mu);
		  su3 ABC,ABCF;
		  unsafe_su3_prod_su3(ABC,conf[par][A.nastyConvert()][nu.nastyConvert()],conf[(1-par).nastyConvert()][B.nastyConvert()][mu.nastyConvert()]);
		  unsafe_su3_prod_su3_dag(ABCF,ABC,conf[(1-par)][F.nastyConvert()][nu.nastyConvert()]);
		  
		  //taking the trace we summ to plaq_summ (only if nu>mu)
		  if(nu>mu)
		    point_shapes[ivol.nastyConvert()][RE]+=real_part_of_trace_su3_prod_su3_dag(ABCF,conf[par][A.nastyConvert()][mu.nastyConvert()]);
		  
		  //compute backward staple starting from A
		  su3 ADE,ADEF;
		  unsafe_su3_dag_prod_su3(ADE,conf[(1-par).nastyConvert()][D.nastyConvert()][nu.nastyConvert()],conf[(1-par)][D.nastyConvert()][mu.nastyConvert()]);
		  unsafe_su3_prod_su3(ADEF,ADE,conf[par][E.nastyConvert()][nu.nastyConvert()]);
		  
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
    
    FOR_ALL_DIRS(mu) //link dir
      FOR_ALL_DIRS(nu) //staple dir
	if(nu!=mu)
	  {
	    NISSA_PARALLEL_LOOP(A,0,locVol)
	      {
		//compute forward staple starting from A
		const LocLxSite& B=loclxNeighup(A,nu),D=loclxNeighdw(A,nu);
		const LocLxSite& E=loclxNeighup(D,mu),F=loclxNeighup(A,mu);
		su3 ABC,ABCF;
		unsafe_su3_prod_su3(ABC,conf[A.nastyConvert()][nu.nastyConvert()],conf[B.nastyConvert()][mu.nastyConvert()]);
		unsafe_su3_prod_su3_dag(ABCF,ABC,conf[F.nastyConvert()][nu.nastyConvert()]);
		
		//taking the trace we summ to plaq_summ (only if nu>mu)
		if(nu>mu) point_shapes[A.nastyConvert()][RE]+=real_part_of_trace_su3_prod_su3_dag(ABCF,conf[A.nastyConvert()][mu.nastyConvert()]);
		
		//compute backward staple starting from A
		su3 ADE,ADEF;
		unsafe_su3_dag_prod_su3(ADE,conf[D.nastyConvert()][nu.nastyConvert()],conf[D.nastyConvert()][mu.nastyConvert()]);
		unsafe_su3_prod_su3(ADEF,ADE,conf[E.nastyConvert()][nu.nastyConvert()]);
		
		//taking the trace we summ to rect_summ
		point_shapes[A.nastyConvert()][IM]+=real_part_of_trace_su3_prod_su3_dag(ABCF,ADEF);
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
    complex *loc_shapes=nissa_malloc("loc_shapes",glbTimeSize.nastyConvert(),complex);
    vector_reset(loc_shapes);
    
    //loop over time
    NISSA_PARALLEL_LOOP(loc_t,0,locTimeSize)
      for(LocLxSite ivol=loc_t*locSpatVol;ivol<(loc_t+1)*locSpatVol;ivol++)
	complex_summassign(loc_shapes[glbCoordOfLoclx(ivol,tDir).nastyConvert()],point_shapes[ivol.nastyConvert()]);
    NISSA_PARALLEL_LOOP_END;
    nissa_free(point_shapes);
    
    //reduce (passing throug additional var because of external unkwnon env)
    complex *coll_shapes=nissa_malloc("coll_shapes",glbTimeSize.nastyConvert(),complex);
    if(IS_MASTER_THREAD) MPI_Reduce(loc_shapes,coll_shapes,2*glbTimeSize(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    nissa_free(loc_shapes);
    
    //normalize
    FOR_ALL_GLB_TIMES(t)
      {
	glb_shapes[2*t.nastyConvert()+0]=coll_shapes[t.nastyConvert()][RE]/(18.0*glbVol()/glbTimeSize());
	glb_shapes[2*t.nastyConvert()+1]=coll_shapes[t.nastyConvert()][IM]/(36.0*glbVol()/glbTimeSize());
      }
    nissa_free(coll_shapes);
  }
}
