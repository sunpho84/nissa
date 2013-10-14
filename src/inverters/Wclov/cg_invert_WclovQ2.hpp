#ifndef _CG_INVERT_WCLOVQ2_H
#define _CG_INVERT_WCLOVQ2_H

namespace nissa
{
  void inv_WclovQ2_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,int niter,double residue,spincolor *source);
}

#endif
