#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/vectors.hpp"

#include "dirac_operators/WclovQ/dirac_operator_WclovQ.hpp"
#include "cg_invert_WclovQ2.hpp"

namespace nissa
{
  void inv_WclovQ_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,int niter,double residue,spincolor *source)
  {
    inv_WclovQ2_cg(sol,NULL,conf,kappa,csw,Pmunu,niter,residue,source);
    spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
    vector_copy(temp,sol);
    apply_WclovQ(sol,conf,kappa,csw,Pmunu,temp);
    nissa_free(temp);
  }
}
