#ifndef _CG_INVERT_TMCLOVQ2_64_HPP
#define _CG_INVERT_TMCLOVQ2_64_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmclovQ2_cg_64(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double mu,int niter,double residue,spincolor *source);
}

#endif
