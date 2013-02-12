#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../new_types/new_types_definitions.h"
#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/debug.h"
#include "../../base/vectors.h"

#include "../dirac_operator_tmQ/dirac_operator_tmQ_128.h"

//WORKARAOUND
int allocated_128=0;

void apply_tmQ2_RL_128(spincolor_128 *out,quad_su3 *conf,double kappa,spincolor_128 *temp,int RL,double mu,spincolor_128 *in)
{
#pragma omp single
  if(temp==NULL)
    {
      temp=nissa_malloc("tempQ",loc_vol+bord_vol,spincolor_128);
      allocated_128=1;
    }

#pragma omp single
  if(RL==1) crash("left not implmented");

  apply_tmQ_128(temp,conf,kappa,+mu,in);
  apply_tmQ_128(out,conf,kappa,-mu,temp);

#pragma omp single
  if(allocated_128==1)
    {
      nissa_free(temp);
      temp=NULL;
    }
}
