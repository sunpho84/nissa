#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../base/communicate.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../routines/openmp.h"

#include "dirac_operator_tmQ_portable.cpp"

#include "../dirac_operator_tmQ_left/dirac_operator_tmQ_left.h"

//wrapper - to be moved elsewhere
void apply_tmQ_RL(spincolor *out,quad_su3 *conf,double kappa,double mu,int RL,spincolor *in)
{
  if(RL==0) apply_tmQ(out,conf,kappa,mu,in);
  else apply_tmQ_left(out,conf,kappa,mu,in);

  set_borders_invalid(out);
}
