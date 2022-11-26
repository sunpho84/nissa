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
  void Symanzik_force_lx_conf(LxField<quad_su3>& out,
			      const LxField<quad_su3>& conf,
			      const double& beta,
			      const double& C1)
  {
    verbosity_lv2_master_printf("Computing Symanzik force\n");
    
    //coefficient of rectangles and squares, including beta
    const double C0=get_C0(C1);
    const double w1=-C1*beta/NCOL,w0=-C0*beta/NCOL;
    
    
    //compute squared pieces
    LxField<squared_staples_t> squared_staples("squared_staples",WITH_HALO);
    compute_squared_staples_lx_conf(squared_staples,conf);
    
    //compute rectangular pieces
    LxField<rectangular_staples_t> rectangular_staples("rectangular_staples",WITH_HALO);
    compute_rectangular_staples_lx_conf(rectangular_staples,conf,squared_staples);
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  //summ the six terms of squares
	  su3_summ(out[ivol][mu],squared_staples[ivol][mu][0],squared_staples[ivol][mu][1]);
	  for(int iterm=2;iterm<6;iterm++) su3_summassign(out[ivol][mu],squared_staples[ivol][mu][iterm]);
	  safe_su3_hermitian_prod_double(out[ivol][mu],out[ivol][mu],w0);
	  
	  //summ the six terms of rectangles
	  su3 temp;
	  su3_summ(temp,rectangular_staples[ivol][mu][0],rectangular_staples[ivol][mu][1]);
	  for(int iterm=2;iterm<6;iterm++) su3_summassign(temp,rectangular_staples[ivol][mu][iterm]);
	  su3_summ_the_hermitian_prod_double(out[ivol][mu],temp,w1);
	}
    NISSA_PARALLEL_LOOP_END;
    
    out.invalidateHalo();
  }
}
