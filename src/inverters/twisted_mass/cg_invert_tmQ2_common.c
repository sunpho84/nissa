#pragma once

void inv_tmQ2_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,int rniter,double residue,spincolor *source)
{
  inv_tmQ2_cg_RL(sol,guess,conf,kappa,m,niter,rniter,residue,0,source);
}

void inv_tmQ2_cg_left(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,int rniter,double residue,spincolor *source)
{
  inv_tmQ2_cg_RL(sol,guess,conf,kappa,m,niter,rniter,residue,1,source);
}
