#ifndef _CG_INVERT_TMCLOVQ2_64_BGQ_H
#define _CG_INVERT_TMCLOVQ2_64_BGQ_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void inv_tmclovQ2_cg_64_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 *conf,double kappa,bi_opt_as2t_su3 *Cl,double mu,int niter,double residue,bi_spincolor *source);
}

#endif
