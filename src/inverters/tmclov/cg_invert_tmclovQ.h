#ifndef _CG_INVERT_TMCLOVQ_H
#define _CG_INVERT_TMCLOVQ_H

void inv_tmclovQ_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double mu,int niter,int rniter,double residue,spincolor *source);

#endif
