#ifndef _CG_INVERT_TMCLOVQ2_64_HPP
#define _CG_INVERT_TMCLOVQ2_64_HPP

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"

#include "cg_64_invert_tmclovQ2_portable.hpp"

namespace nissa
{
  inline void inv_tmclovQ2_cg_64(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,clover_term_t *Cl,double mu,int niter,double residue,spincolor *source)
  {
    inv_tmclovQ2_cg_64_portable(sol,guess,conf,kappa,Cl,mu,niter,residue,source);
  }
}

#endif
