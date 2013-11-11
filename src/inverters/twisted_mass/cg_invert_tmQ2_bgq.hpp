#ifndef _CG_INVERT_TMQ2_BGQ_H
#define _CG_INVERT_TMQ2_BGQ_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void inv_tmQ2_RL_cg_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 *conf,double kappa,int RL,double m,int niter,double residue,bi_spincolor *source);
  void inv_tmQ2_cg_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 *conf,double kappa,double m,int niter,double residue,bi_spincolor *source);
  void inv_tmQ2_cg_left_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 *conf,double kappa,double m,int niter,double residue,bi_spincolor *source);
}

#endif
