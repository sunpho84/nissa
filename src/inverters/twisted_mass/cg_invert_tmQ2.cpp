#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/float_128.hpp"
#include "cg_64_invert_tmQ2.hpp"
#include "cg_128_invert_tmQ2.hpp"

namespace nissa
{
  //switch 64 and 128
  void inv_tmQ2_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,double residue,spincolor *source)
  {
      crash("reimplement");//not linking
      //if(use_128_bit_precision)
    //inv_tmQ2_cg_128(sol,guess,conf,kappa,m,niter,residue,source);
      //else
      //inv_tmQ2_cg_64(sol,guess,conf,kappa,m,niter,residue,source);
  }
}
