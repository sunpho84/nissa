#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "operations/su3_paths/rectangular_staples.hpp"
#include "operations/su3_paths/squared_staples.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
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
    
    //compute rectangular pieces
    rectangular_staples_t *rectangular_staples=nissa_malloc("rectangular_staples",loc_vol+bord_vol,rectangular_staples_t);
    compute_rectangular_staples_lx_conf(rectangular_staples,conf,squared_staples);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<4;mu++)
	{
	  //summ the six terms of squares
	  su3_summ(out[ivol][mu],squared_staples[ivol][mu][0],squared_staples[ivol][mu][1]);
	  for(int iterm=2;iterm<6;iterm++) su3_summassign(out[ivol][mu],squared_staples[ivol][mu][iterm]);
	  safe_su3_hermitian_prod_double(out[ivol][mu],out[ivol][mu],c0);

	  //summ the six terms of rectangles
	  su3 temp;
	  su3_summ(temp,rectangular_staples[ivol][mu][0],rectangular_staples[ivol][mu][1]);
	  for(int iterm=2;iterm<6;iterm++) su3_summassign(temp,rectangular_staples[ivol][mu][iterm]);
	  su3_summ_the_hermitian_prod_double(out[ivol][mu],temp,c1);
	}
    
    nissa_free(squared_staples);
    nissa_free(rectangular_staples);
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
}
