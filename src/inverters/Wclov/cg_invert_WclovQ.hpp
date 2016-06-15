#ifndef _CG_INVERT_WCLOVQ_HPP
#define _CG_INVERT_WCLOVQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_WclovQ_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,clover_term_t *Cl,int niter,double residue,spincolor *source);
}

#endif
