#ifndef _CG_128_INVERT_TMQ2_H
#define _CG_128_INVERT_TMQ2_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void inv_tmQ2_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double mass,int niter,int rniter,double external_solver_residue,spincolor *external_source);
  void inv_tmQ2_m2_RL_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,int RL,double m2,int niter,int rniter,double external_solver_residue,spincolor *external_source);
  void inv_tmQ2_RL_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,int RL,double mass,int niter,int rniter,double external_solver_residue,spincolor *external_source);
}

#endif
