#pragma once

#ifdef BGP

#include "dirac_operator_tmQ_bgp.cpp"

#else

#include "dirac_operator_tmQ_portable.cpp"

#endif

#include "dirac_operator_tmQ_128_portable.cpp"

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

