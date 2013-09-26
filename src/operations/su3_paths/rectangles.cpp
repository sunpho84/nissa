#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //compute plaquettes and rectangles
  THREADABLE_FUNCTION_2ARG(global_plaquette_and_rectangles_eo_conf, double*,glb_shapes, quad_su3**,conf)
  {
    GET_THREAD_ID();
    
    communicate_eo_quad_su3_edges(conf);
    
    //summ squares and rectangles separately
    complex *point_shapes=nissa_malloc("point_shapes",loc_vol,complex);
    vector_reset(point_shapes);
    
    for(int par=0;par<2;par++)
      for(int mu=0;mu<4;mu++) //link dir
	for(int nu=0;nu<4;nu++) //staple dir
	  if(nu!=mu)
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
    THREAD_BARRIER();
    
    //reduce and free
    complex coll_shapes;
    complex_vector_glb_collapse(coll_shapes,point_shapes,loc_vol);
    nissa_free(point_shapes);
    
    //normalize (passing throug additional var because of external unkwnon env)
    glb_shapes[RE]=coll_shapes[RE]/(18*glb_vol);
    glb_shapes[IM]=coll_shapes[IM]/(36*glb_vol);
  }}
}
