#pragma once

void inv_Q2_cg(spincolor *sol,spincolor *source,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,int rniter,double residue)
{
  inv_Q2_cg_RL(sol,source,guess,conf,kappa,m,niter,rniter,residue,0);
}

void inv_Q2_cg_left(spincolor *sol,spincolor *source,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,int rniter,double residue)
{
  inv_Q2_cg_RL(sol,source,guess,conf,kappa,m,niter,rniter,residue,1);
}
