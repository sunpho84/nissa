#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/float_128.hpp"
#include "cg_64_invert_tmclovQ2.hpp"
#include "cg_128_invert_tmclovQ2.hpp"

#include "geometry/geometry_lx.hpp"

namespace nissa
{
  //switch 64 and 128
  void inv_tmclovQ2_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,clover_term_t *Cl,double m,int niter,double residue,spincolor *source)
  {
    if(use_128_bit_precision) inv_tmclovQ2_cg_128(sol,guess,conf,kappa,Cl,m,niter,residue,source);
    else inv_tmclovQ2_cg_64(sol,guess,conf,kappa,Cl,m,niter,residue,source);
  }
}
