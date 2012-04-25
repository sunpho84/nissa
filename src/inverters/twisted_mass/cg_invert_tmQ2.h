#ifndef _CG_INVERT_TMQ2_H
#define _CG_INVERT_TMQ2_H

void inv_tmQ2_RL_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,int RL,double m,int niter,int rniter,double residue,spincolor *source);
void inv_tmQ2_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,int rniter,double residue,spincolor *source);
void inv_tmQ2_cg_left(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,int rniter,double residue,spincolor *source);

#endif
