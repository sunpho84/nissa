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
