#pragma once

#include "dirac_operator_tmQ/dirac_operator_tmQ.c"
#include "dirac_operator_tmQ_left/dirac_operator_tmQ_left.c"
#include "dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr.c"

#include "dirac_operator_stD/dirac_operator_stD.c"

///////////////////////////////////////////// TWISTED MASS DEGENERATE LIGHT DOUBLET ///////////////////////////////////////

//wrapper
void apply_tmQ_RL(spincolor *out,quad_su3 *conf,double kappa,double mu,int RL,spincolor *in)
{
  if(RL==0) apply_tmQ(out,conf,kappa,mu,in);
  else apply_tmQ_left(out,conf,kappa,mu,in);

  set_borders_invalid(out);
}

//Apply the Q+ and Q- operator to a spincolor,so that we have Q-^-1 (r==0) and Q+^-1 (r==1) as output
void reconstruct_tm_doublet(spincolor *outminus,spincolor *outplus,quad_su3 *conf,double kappac,double mu,spincolor *in)
{
  apply_tmQ(outminus,conf,kappac,mu,in);
  nissa_loc_vol_loop(ivol)
    unsafe_spincolor_summ_with_ifactor(outplus[ivol],outminus[ivol],in[ivol],-2*mu);
  
  set_borders_invalid(outminus);
  set_borders_invalid(outplus);
}

//Apply the Q+Q- operator to a spincolor
void apply_tmQ2_RL(spincolor *out,quad_su3 *conf,double kappa,double mu,spincolor *temp,int RL,spincolor *in)
{
  int all=0;

  if(temp==NULL)
    {
      temp=nissa_malloc("tempQ",loc_vol+bord_vol,spincolor);
      all=1;
    }

  if(RL==0) apply_tmQ(temp,conf,kappa,+mu,in);
  else apply_tmQ_left(temp,conf,kappa,+mu,in);
  
  communicate_lx_spincolor_borders(temp);
  if(RL==0) apply_tmQ(out,conf,kappa,-mu,temp);
  else apply_tmQ_left(out,conf,kappa,-mu,temp);

  if(all==1)
    {
      nissa_free(temp);
      temp=NULL;
    }
  
  set_borders_invalid(out);
}

void apply_tmQ2(spincolor *out,quad_su3 *conf,double kappa,double mu,spincolor *temp,spincolor *in)
{apply_tmQ2_RL(out,conf,kappa,mu,temp,0,in);}

void apply_tmQ2_left(spincolor *out,quad_su3 *conf,double kappa,double mu,spincolor *temp,spincolor *in)
{apply_tmQ2_RL(out,conf,kappa,mu,temp,1,in);}
