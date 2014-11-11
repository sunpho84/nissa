#ifndef _CG_64_INVERT_TMDEOIMPR_BGQ_HPP
#define _CG_64_INVERT_TMDEOIMPR_BGQ_HPP

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void inv_tmDkern_eoprec_square_eos_cg_64_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 **bi_eo_conf,double kappa,double mu,int niter,double residue,bi_spincolor *source);
}

#endif
