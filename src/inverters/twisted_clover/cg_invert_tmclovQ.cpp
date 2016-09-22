#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"

#include "dirac_operators/tmclovQ/dirac_operator_tmclovQ.hpp"
#include "cg_invert_tmclovQ2.hpp"

namespace nissa
{
  void inv_tmclovQ_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,clover_term_t *Cl,double mu,int niter,double residue,spincolor *source)
  {
    inv_tmclovQ2_cg(sol,NULL,conf,kappa,Cl,mu,niter,residue,source);
    spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
    
    //remove the "wrong r"
    vector_copy(temp,sol);
    apply_tmclovQ(sol,conf,kappa,Cl,-mu,temp);
    
    nissa_free(temp);
  }
}
