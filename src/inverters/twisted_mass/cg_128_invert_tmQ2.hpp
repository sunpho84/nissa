#ifndef _CG_128_INVERT_TMQ2_HPP
#define _CG_128_INVERT_TMQ2_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmQ2_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double mass,int niter,double external_solver_residue,spincolor *external_source);
  void inv_tmQ2_m2_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double m2,int niter,double external_solver_residue,spincolor *external_source);
  void inv_tmQ2_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double mass,int niter,double external_solver_residue,spincolor *external_source);
}

#endif
