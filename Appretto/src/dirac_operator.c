#pragma once

//Include the dirac operator: this depend on the machine
#ifdef BGP

#include "dirac_operator_bgp.c"

#elif defined SSE

#else

#include "dirac_operator_portable.c"

#endif

//Apply the Q+ and Q- operator to a spincolor
void reconstruct_doublet(spincolor *outplus,spincolor *outminus,spincolor *in,quad_su3 *conf,double kappac,double mu)
{
  apply_Q(outplus,in,conf,kappac,mu);
  for(int X=0;X<loc_vol;X++)
    unsafe_spincolor_summ_with_ifactor(outminus[X],outplus[X],in[X],-2*mu);
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
