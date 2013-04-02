#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/communicate.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../new_types/new_types_definitions.h"
#include "../../routines/ios.h"
#include "../../routines/openmp.h"

//compute the tree level Symanzik action
THREADABLE_FUNCTION_4ARG(tree_level_Symanzik_action, double*,action, quad_su3**,conf, double,beta, int,stagphases_present)
{
  GET_THREAD_ID();
  
  verbosity_lv1_master_printf("Computing tree level Symanzik action\n");
  
  //coefficient of rectangles and squares
  double b1=-1.0/12,b0=1-8*b1;
  b0=0;
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
  
  //reduce and free
  complex glb_shapes;
  complex_vector_glb_collapse(glb_shapes,point_shapes,loc_vol);
  nissa_free(point_shapes);

  if(stagphases_present) glb_shapes[RE]*=-1; //stag phases add (-1)^area
  
  //compute the total action
  (*action)=(b0*(18*glb_vol-glb_shapes[RE])+b1*(36*glb_vol-glb_shapes[IM]))*beta/3;
}}
