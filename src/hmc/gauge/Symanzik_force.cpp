#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "operations/su3_paths/rectangular_staples.hpp"
#include "operations/su3_paths/squared_staples.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

#include "Symanzik_action.hpp"

namespace nissa
{
  void Symanzik_force_lx_conf(quad_su3* out,quad_su3* conf,double beta,double C1)
  {
    verbosity_lv2_master_printf("Computing Symanzik force\n");
    
    //coefficient of rectangles and squares, including beta
    double C0=get_C0(C1);
    double w1=-C1*beta/NCOL,w0=-C0*beta/NCOL;
    
    
    //compute squared pieces
    squared_staples_t *squared_staples=nissa_malloc("squared_staples",locVolWithBord.nastyConvert(),squared_staples_t);
    compute_squared_staples_lx_conf(squared_staples,conf);
    
    //compute rectangular pieces
    rectangular_staples_t *rectangular_staples=nissa_malloc("rectangular_staples",locVolWithBord.nastyConvert(),rectangular_staples_t);
    compute_rectangular_staples_lx_conf(rectangular_staples,conf,squared_staples);
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      FOR_ALL_DIRS(mu)
	{
	  //summ the six terms of squares
	  su3_summ(out[ivol.nastyConvert()][mu.nastyConvert()],squared_staples[ivol.nastyConvert()][mu.nastyConvert()][0],squared_staples[ivol.nastyConvert()][mu.nastyConvert()][1]);
	  for(int iterm=2;iterm<6;iterm++) su3_summassign(out[ivol.nastyConvert()][mu.nastyConvert()],squared_staples[ivol.nastyConvert()][mu.nastyConvert()][iterm]);
	  safe_su3_hermitian_prod_double(out[ivol.nastyConvert()][mu.nastyConvert()],out[ivol.nastyConvert()][mu.nastyConvert()],w0);
	  
	  //summ the six terms of rectangles
	  su3 temp;
	  su3_summ(temp,rectangular_staples[ivol.nastyConvert()][mu.nastyConvert()][0],rectangular_staples[ivol.nastyConvert()][mu.nastyConvert()][1]);
	  for(int iterm=2;iterm<6;iterm++) su3_summassign(temp,rectangular_staples[ivol.nastyConvert()][mu.nastyConvert()][iterm]);
	  su3_summ_the_hermitian_prod_double(out[ivol.nastyConvert()][mu.nastyConvert()],temp,w1);
	}
    NISSA_PARALLEL_LOOP_END;
    
    nissa_free(squared_staples);
    nissa_free(rectangular_staples);
    
    set_borders_invalid(out);
  }
}
