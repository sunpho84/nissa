#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/vectors.hpp"

#include "dirac_operators/tmclovQ/dirac_operator_tmclovQ.hpp"
#include "cg_invert_tmclovQ2.hpp"

namespace nissa
{
  void inv_tmclovQ_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double mu,int niter,double residue,spincolor *source)
  {
    inv_tmclovQ2_cg(sol,NULL,conf,kappa,csw,Pmunu,mu,niter,residue,source);
    spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
    
    vector_copy(temp,sol);
    apply_tmclovQ(sol,conf,kappa,csw,Pmunu,mu,temp);
    
    nissa_free(temp);
  }
}
