#pragma once

//Include the dirac operator: this depend on the machine
#ifdef BGP

#include "dirac_operator_bgp.c"

#elif defined SSE
//possibly to be added
#else

#include "dirac_operator_portable.c"

#endif

//Apply the Q+ and Q- operator to a spincolor,so that we have Q-^-1 (r==0) and Q+^-1 (r==1) as output
void reconstruct_doublet(spincolor *outminus,spincolor *outplus,spincolor *in,quad_su3 *conf,double kappac,double mu)
{
  apply_Q(outminus,in,conf,kappac,mu);
  for(int X=0;X<loc_vol;X++)
    unsafe_spincolor_summ_with_ifactor(outplus[X],outminus[X],in[X],-2*mu);
}

//Apply the Q+Q- operator to a spincolor
void apply_Q2(spincolor *out,spincolor *in,quad_su3 *conf,double kappa,double mu,spincolor *temp)
{
  int all=0;

  if(temp==NULL)
    {
      temp=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
      all=1;
    }

  apply_Q(temp,in,conf,kappa,+mu);
  communicate_lx_spincolor_borders(temp);
  apply_Q(out,temp,conf,kappa,-mu);

  if(all==1)
    {
      free(temp);
      temp=NULL;
    }
}
