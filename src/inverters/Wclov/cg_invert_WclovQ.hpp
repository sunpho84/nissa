#ifndef _CG_INVERT_WCLOVQ_H
#define _CG_INVERT_WCLOVQ_H

namespace nissa
{
  void inv_WclovQ_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,int niter,int rniter,double residue,spincolor *source);
}

#endif
