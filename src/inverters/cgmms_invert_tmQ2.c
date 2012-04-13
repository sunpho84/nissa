#pragma once


void inv_tmQ2_RL_cgmms(spincolor **sol,quad_su3 *conf,double kappa,int RL,double *m,int nmass,int niter_max,double req_res,spincolor *source)
{
  double m2[nmass];
  for(int imass=0;imass<nmass;imass++) m2[imass]=m[imass]*m[imass];
  inv_tmQ2_RL_cgmm2s(sol,conf,kappa,RL,m2,nmass,niter_max,req_res,source);
}
void inv_tmQ2_cgmms(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double req_res,spincolor *source)
{inv_tmQ2_RL_cgmms(sol,conf,kappa,0,m,nmass,niter_max,req_res,source);}
void inv_tmQ2_cgmms_left(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double req_res,spincolor *source)
{inv_tmQ2_RL_cgmms(sol,conf,kappa,1,m,nmass,niter_max,req_res,source);}
