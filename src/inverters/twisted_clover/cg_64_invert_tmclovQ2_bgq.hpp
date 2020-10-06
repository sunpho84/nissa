#ifndef _CG_INVERT_TMCLOVQ2_64_BGQ_HPP
#define _CG_INVERT_TMCLOVQ2_64_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmclovQ2_cg_64_bgq(vir_spincolor *sol,vir_spincolor *guess,vir_oct_su3 *conf,double kappa,vir_clover_term_t *Cl,double mu,int niter,double residue,vir_spincolor *source);
}

#endif
