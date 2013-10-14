#ifndef _CG_INVERT_TMCLOVQ2_H
#define _CG_INVERT_TMCLOVQ2_H

namespace nissa
{
  void inv_tmclovQ2_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double mu,int niter,double residue,spincolor *source);
}

#endif
