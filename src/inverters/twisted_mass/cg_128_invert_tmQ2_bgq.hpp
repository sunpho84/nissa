#ifndef _CG_128_INVERT_TMQ2_BGQ_HPP
#define _CG_128_INVERT_TMQ2_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmQ2_cg_128_bgq(vir_spincolor *sol,vir_spincolor *guess,vir_oct_su3 *conf,double kappa,double mass,int niter,double external_solver_residue,vir_spincolor *external_source);
  void inv_tmQ2_m2_RL_cg_128_bgq(vir_spincolor *sol,vir_spincolor *guess,vir_oct_su3 *conf,double kappa,int RL,double m2,int niter,double external_solver_residue,vir_spincolor *external_source);
  void inv_tmQ2_RL_cg_128_bgq(vir_spincolor *sol,vir_spincolor *guess,vir_oct_su3 *conf,double kappa,int RL,double mass,int niter,double external_solver_residue,vir_spincolor *external_source);
}

#endif
