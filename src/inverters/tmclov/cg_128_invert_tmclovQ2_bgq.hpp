#ifndef _CG_128_INVERT_TMCLOVQ2_BGQ_H
#define _CG_128_INVERT_TMCLOVQ2_BGQ_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void inv_tmclovQ2_cg_128_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 *conf,double kappa,bi_opt_as2t_su3 *Cl,double mass,int niter,double external_solver_residue,bi_spincolor *external_source);
  void inv_tmclovQ2_m2_cg_128_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 *conf,double kappa,bi_opt_as2t_su3 *Cl,double m2,int niter,double external_solver_residue,bi_spincolor *external_source);
}

#endif
