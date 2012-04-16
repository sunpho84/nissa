#pragma once

void inv_tmQ2_m2_RL_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,int RL,double m2,int niter,double external_solver_residue,spincolor *external_source)
{inv_tmQ2_RL_cg_128(sol,guess,conf,kappa,RL,sqrt(m2),niter,external_solver_residue,external_source);}
