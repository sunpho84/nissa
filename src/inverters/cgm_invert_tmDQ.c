#pragma once

////////////////////////////////////// TWISTED MASS LIGHT DEGENERATE INVERTERS ///////////////////////////////

void inv_tmDQ_RL_cgm(spincolor **sol,quad_su3 *conf,double kappa,int RL,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
{
  //put the g5
  nissa_loc_vol_loop(ivol) for(int id1=2;id1<4;id1++) for(int ic1=0;ic1<3;ic1++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic1][ri]*=-1;
  set_borders_invalid(source);
  inv_tmQ2_RL_cgm(sol,conf,kappa,RL,m,nmass,niter_max,req_res,source);
  nissa_loc_vol_loop(ivol) for(int id1=2;id1<4;id1++) for(int ic1=0;ic1<3;ic1++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic1][ri]*=-1;
  set_borders_invalid(source);
}

void inv_tmDQ_cgm(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
{inv_tmDQ_RL_cgm(sol,conf,kappa,0,m,nmass,niter_max,req_res,source);}
void inv_tmDQ_cgm_left(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
{inv_tmDQ_RL_cgm(sol,conf,kappa,1,m,nmass,niter_max,req_res,source);}

