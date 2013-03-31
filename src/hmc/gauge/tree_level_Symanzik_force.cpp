#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <string.h>

#include "../../base/communicate.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../geometry/geometry_eo.h"
#include "../../geometry/geometry_mix.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../operations/su3_paths/plaquette.h"
#include "../../operations/su3_paths/arbitrary.h"
#include "../../routines/ios.h"

THREADABLE_FUNCTION_3ARG(tree_level_Symanzik_force_lx_conf, quad_su3*,out, quad_su3*,conf, double,beta)
{
  verbosity_lv1_master_printf("Computing tree level Symanzik force\n");
  
  //coefficient of rectangles and squares, including beta
  double b1=-1.0/12,b0=1-8*b1;
  double c1=-b1*beta/3,c0=b0*beta/3; //the stag phases add (-1)^area
  
  GET_THREAD_ID();
  
  //compute squared pieces
  squared_staples_t *squared_staples=nissa_malloc("squared_staples",loc_vol+bord_vol,squared_staples_t);
  compute_squared_staples_lx_conf(squared_staples,conf);

  //compute the squared staples
  compute_summed_squared_staples_lx_conf(out,conf);

  //take hermitian*r
  double r=beta/3;
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<4;mu++)
      safe_su3_hermitian_prod_double(out[ivol][mu],out[ivol][mu],r);
}}
