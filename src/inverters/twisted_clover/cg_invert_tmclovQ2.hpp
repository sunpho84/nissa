#ifndef _CG_INVERT_TMCLOVQ2_HPP
#define _CG_INVERT_TMCLOVQ2_HPP

namespace nissa
{
  void inv_tmclovQ2_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,clover_term_t *Cl,double mu,int niter,double residue,spincolor *source);
}

#endif
