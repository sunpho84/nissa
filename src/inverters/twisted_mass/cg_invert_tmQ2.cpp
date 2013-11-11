#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "cg_64_invert_tmQ2.hpp"
#include "cg_128_invert_tmQ2.hpp"

#include "new_types/new_types_definitions.hpp"
#include "base/global_variables.hpp"

namespace nissa
{
  //switch 64 and 128
  void inv_tmQ2_RL_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,int RL,double m,int niter,double residue,spincolor *source)
  {
    if(use_128_bit_precision) inv_tmQ2_RL_cg_128(sol,guess,conf,kappa,RL,m,niter,residue,source);
    else inv_tmQ2_RL_cg_64(sol,guess,conf,kappa,RL,m,niter,residue,source);
  }
}
