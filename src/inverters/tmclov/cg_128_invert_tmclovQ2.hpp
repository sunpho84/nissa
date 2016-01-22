#ifndef _CG_128_INVERT_TMCLOVQ2_HPP
#define _CG_128_INVERT_TMCLOVQ2_HPP

namespace nissa
{
  void inv_tmclovQ2_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double mass,int niter,double external_solver_residue,spincolor *external_source);
  void inv_tmclovQ2_m2_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double m2,int niter,double external_solver_residue,spincolor *external_source);
}

#endif
