#ifndef _CG_64_INVERT_TMDEOIMPR_BGQ_HPP
#define _CG_64_INVERT_TMDEOIMPR_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmDkern_eoprec_square_eos_cg_64_bgq(vir_spincolor *sol,vir_spincolor *guess,vir_oct_su3 **vir_eo_conf,double kappa,double mu,int niter,double residue,vir_spincolor *source);
}

#endif
