#pragma once

//Include the dirac operator: this depend on the machine
#ifdef BGP

#include "dirac_operator_tmQ_bgp.c"
#include "dirac_operator_tmQ_left_bgp.c"

#else

#include "dirac_operator_stD_portable.c"

#include "dirac_operator_tmQ_portable.c"
#include "dirac_operator_tmQ_left_portable.c"

#endif

///////////////////////////////////////////// TWISTED MASS DEGENERATE LIGHT DOUBLET ///////////////////////////////////////

//wrapper
void apply_tmQ_RL(spincolor *out,spincolor *in,quad_su3 *conf,double kappa,double mu,int RL)
{
  if(RL==0) apply_tmQ(out,in,conf,kappa,mu);
  else apply_tmQ_left(out,in,conf,kappa,mu);
}

//Apply the Q+ and Q- operator to a spincolor,so that we have Q-^-1 (r==0) and Q+^-1 (r==1) as output
void reconstruct_tm_doublet(spincolor *outminus,spincolor *outplus,spincolor *in,quad_su3 *conf,double kappac,double mu)
{
  apply_tmQ(outminus,in,conf,kappac,mu);
  for(int X=0;X<loc_vol;X++)
    unsafe_spincolor_summ_with_ifactor(outplus[X],outminus[X],in[X],-2*mu);
}

//Apply the Q+Q- operator to a spincolor
void apply_tmQ2_RL(spincolor *out,spincolor *in,quad_su3 *conf,double kappa,double mu,spincolor *temp,int RL)
{
  int all=0;

  if(temp==NULL)
    {
      temp=nissa_malloc("tempQ",loc_vol+loc_bord,spincolor);
      all=1;
    }

  if(RL==0) apply_tmQ(temp,in,conf,kappa,+mu);
  else apply_tmQ_left(temp,in,conf,kappa,+mu);
  
  communicate_lx_spincolor_borders(temp);
  if(RL==0) apply_tmQ(out,temp,conf,kappa,-mu);
  else apply_tmQ_left(out,temp,conf,kappa,-mu);

  if(all==1)
    {
      nissa_free(temp);
      temp=NULL;
    }
}

void apply_tmQ2(spincolor *out,spincolor *in,quad_su3 *conf,double kappa,double mu,spincolor *temp)
{apply_tmQ2_RL(out,in,conf,kappa,mu,temp,0);}

void apply_tmQ2_left(spincolor *out,spincolor *in,quad_su3 *conf,double kappa,double mu,spincolor *temp)
{apply_tmQ2_RL(out,in,conf,kappa,mu,temp,1);}
