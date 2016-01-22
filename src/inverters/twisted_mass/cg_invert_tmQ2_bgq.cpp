#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "cg_64_invert_tmQ2_bgq.hpp"
#include "cg_128_invert_tmQ2_bgq.hpp"

namespace nissa
{
  //switch 64 and 128
  void inv_tmQ2_RL_cg_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 *conf,double kappa,int RL,double m,int niter,double residue,bi_spincolor *source)
  {
    if(use_128_bit_precision) inv_tmQ2_RL_cg_128_bgq(sol,guess,conf,kappa,RL,m,niter,residue,source);
    else inv_tmQ2_RL_cg_64_bgq(sol,guess,conf,kappa,RL,m,niter,residue,source);
  }
}
