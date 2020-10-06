#ifndef _CG_INVERT_TMQ2_BGQ_HPP
#define _CG_INVERT_TMQ2_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmQ2_RL_cg_64_bgq(vir_spincolor *sol,vir_spincolor *guess,vir_oct_su3 *conf,double kappa,int RL,double m,int niter,double residue,vir_spincolor *source);
  void inv_tmQ2_cg_64_bgq(vir_spincolor *sol,vir_spincolor *guess,vir_oct_su3 *conf,double kappa,double m,int niter,double residue,vir_spincolor *source);
  void inv_tmQ2_cg_left_64_bgq(vir_spincolor *sol,vir_spincolor *guess,vir_oct_su3 *conf,double kappa,double m,int niter,double residue,vir_spincolor *source);
}

#endif
