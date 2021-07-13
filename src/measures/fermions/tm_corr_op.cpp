#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "linalgs/reduce.hpp"

#include "tm_corr_op.hpp"

namespace nissa
{
  void tm_corr_op::ins(spincolor *out,const int igamma,spincolor *in)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	unsafe_dirac_prod_spincolor(out[ivol],base_gamma[igamma],in[ivol]);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
}

