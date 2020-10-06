#pragma once

namespace nissa
{
  void inv_tmQ2_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,double residue,spincolor *source)
  {inv_tmQ2_RL_cg(sol,guess,conf,kappa,0,m,niter,residue,source);}
  
  void inv_tmQ2_cg_left(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,double residue,spincolor *source)
  {inv_tmQ2_RL_cg(sol,guess,conf,kappa,1,m,niter,residue,source);}
}
