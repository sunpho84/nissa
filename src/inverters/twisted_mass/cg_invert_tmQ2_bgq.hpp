#ifndef _CG_INVERT_TMQ2_BGQ_HPP
#define _CG_INVERT_TMQ2_BGQ_HPP

namespace nissa
{
  void inv_tmQ2_RL_cg_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 *conf,double kappa,int RL,double m,int niter,double residue,bi_spincolor *source);
  void inv_tmQ2_cg_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 *conf,double kappa,double m,int niter,double residue,bi_spincolor *source);
  void inv_tmQ2_cg_left_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 *conf,double kappa,double m,int niter,double residue,bi_spincolor *source);
}

#endif
