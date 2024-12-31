#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/field.hpp"
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
  void global_plaquette_and_rectangles_eo_conf(double* glb_shapes,
					       const EoField<quad_su3>& conf)
  {
    conf.updateEdges();
    
    //summ squares and rectangles separately
    LxField<complex> point_shapes("point_shapes");
    point_shapes.reset();
    
    for(int par=0;par<2;par++)
      for(int mu=0;mu<NDIM;mu++) //link dir
	for(int nu=0;nu<NDIM;nu++) //staple dir
	  if(nu!=mu)
	    {
	      PAR(0,
		  locVolh,
		  CAPTURE(par,
			  nu,
			  mu,
			  TO_WRITE(point_shapes),
			  TO_READ(conf)),
		  A,
		  {
		    const int ivol=
		      loclx_of_loceo[par][A];
		    
		    //compute forward staple starting from A
		    const int B=
		      loceo_neighup[par][A][nu];
		    
		    const int D=
		      loceo_neighdw[par][A][nu];
		    
		    const int E=
		      loceo_neighup[!par][D][mu];
		    
		    const int F=
		      loceo_neighup[par][A][mu];
		    
		    su3 ABC,ABCF;
		    unsafe_su3_prod_su3(ABC,conf[par][A][nu],conf[!par][B][mu]);
		    unsafe_su3_prod_su3_dag(ABCF,ABC,conf[!par][F][nu]);
		    
		    //taking the trace we summ to plaq_summ (only if nu>mu)
		    if(nu>mu)
		      point_shapes[ivol][RE]+=real_part_of_trace_su3_prod_su3_dag(ABCF,conf[par][A][mu]);
		    
		    //compute backward staple starting from A
		    su3 ADE,ADEF;
		    unsafe_su3_dag_prod_su3(ADE,conf[!par][D][nu],conf[!par][D][mu]);
		    unsafe_su3_prod_su3(ADEF,ADE,conf[par][E][nu]);
		    
		    //taking the trace we summ to rect_summ
		    point_shapes[ivol][IM]+=real_part_of_trace_su3_prod_su3_dag(ABCF,ADEF);
		  });
	    }
    
    //reduce and free
    complex coll_shapes;
    point_shapes.reduce(coll_shapes);
    
    //normalize (passing throug additional var because of external unkwnon env)
    glb_shapes[RE]=coll_shapes[RE]/(18*glbVol);
    glb_shapes[IM]=coll_shapes[IM]/(36*glbVol);
  }
  
  //compute plaquettes and rectangles
  void point_plaquette_and_rectangles_lx_conf(LxField<complex>& point_shapes,
					      const LxField<quad_su3>& conf)
  {
    //communicate conf and reset point shapes
    conf.updateEdges();
    point_shapes.reset();
    
    for(int mu=0;mu<NDIM;mu++) //link dir
      for(int nu=0;nu<NDIM;nu++) //staple dir
	if(nu!=mu)
	  {
	    PAR(0,locVol,
		CAPTURE(mu,nu,
			TO_WRITE(point_shapes),
			TO_READ(conf)),
		A,
		{
		  const int ivol=A;
		  
		  //compute forward staple starting from A
		  const int B=loclxNeighup[A][nu],D=loclxNeighdw[A][nu];
		  const int E=loclxNeighup[D][mu],F=loclxNeighup[A][mu];
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
		});
	  }
  }
  
  //compute plaquettes and rectangles
  void global_plaquette_and_rectangles_lx_conf(double* glb_shapes,
					       const LxField<quad_su3>& conf)
  {
    //summ squares and rectangles separately
    LxField<complex> point_shapes("point_shapes");
    point_plaquette_and_rectangles_lx_conf(point_shapes,conf);
    
    //reduce and free
    complex coll_shapes;
    glb_reduce(&coll_shapes,point_shapes,locVol);
    
    //normalize (passing throug additional var because of external unkwnon env)
    glb_shapes[RE]=coll_shapes[RE]/(18*glbVol);
    glb_shapes[IM]=coll_shapes[IM]/(36*glbVol);
  }
  
  //compute plaquettes and rectangles
  void global_plaquette_and_rectangles_lx_conf_per_timeslice(double* glb_shapes,
							     const LxField<quad_su3>& conf)
  {
    CRASH("reimplement");
    // //summ squares and rectangles separately
    // LxField<complex> point_shapes("point_shapes");
    // point_plaquette_and_rectangles_lx_conf(point_shapes,conf);
    
    // //reduce
    // complex *loc_shapes=nissa_malloc("loc_shapes",glbSize[0],complex);
    // vector_reset(loc_shapes);
    
    // //loop over time
    // for(int loc_t=0;loc_t<locSize[0];loc_t++)
    //   for(int ivol=loc_t*locSpatVol;ivol<(loc_t+1)*locSpatVol;ivol++)
    // 	complex_summassign(loc_shapes[glbCoordOfLoclx[ivol][0]],point_shapes[ivol]);
    // NISSA_PARALLEL_LOOP_END;
    
    // //reduce (passing throug additional var because of external unkwnon env)
    // complex *coll_shapes=nissa_malloc("coll_shapes",glbSize[0],complex);
    // if(IS_MASTER_THREAD) MPI_Reduce(loc_shapes,coll_shapes,2*glbSize[0],MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    // nissa_free(loc_shapes);
    
    // //normalize
    // for(int t=0;t<glbSize[0];t++)
    //   {
    // 	glb_shapes[2*t+0]=coll_shapes[t][RE]/(18.0*glbVol/glbSize[0]);
    // 	glb_shapes[2*t+1]=coll_shapes[t][IM]/(36.0*glbVol/glbSize[0]);
    //   }
    // nissa_free(coll_shapes);
  }
}
