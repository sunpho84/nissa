#ifndef _CG_INVERT_TMCLOVQ2_64_PORTABLE_HPP
#define _CG_INVERT_TMCLOVQ2_64_PORTABLE_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmclovQ2_cg_64_portable(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,clover_term_t *Cl,double mu,int niter,double residue,spincolor *source);
}

#endif
