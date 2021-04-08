#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"

#include "dirac_operators/WclovQ/dirac_operator_WclovQ.hpp"
#include "cg_invert_WclovQ2.hpp"

namespace nissa
{
  void inv_WclovQ_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,clover_term_t *Cl,int niter,double residue,spincolor *source)
  {
    inv_WclovQ2_cg(sol,NULL,conf,kappa,Cl,niter,residue,source);
    spincolor *temp=nissa_malloc("temp",(locVol+bord_vol).nastyConvert(),spincolor);
    vector_copy(temp,sol);
    apply_WclovQ(sol,conf,kappa,Cl,temp);
    nissa_free(temp);
  }
}
