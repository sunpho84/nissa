#include "../../new_types/new_types_definitions.h"
#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/debug.h"
#include "../../base/vectors.h"

#include "../dirac_operator_tmQ/dirac_operator_tmQ_128.h"


void apply_tmQ2_RL_128(spincolor_128 *out,quad_su3 *conf,double kappa,spincolor_128 *temp,int RL,double mu,spincolor_128 *in)
{
  int all=0;

  if(temp==NULL)
    {
      temp=nissa_malloc("tempQ",loc_vol+bord_vol,spincolor_128);
      all=1;
    }

  if(RL==0) apply_tmQ_128(temp,conf,kappa,+mu,in);
  else crash("left not implmented");
  
  communicate_lx_spincolor_128_borders(temp);
  apply_tmQ_128(out,conf,kappa,-mu,temp);

  if(all==1)
    {
      nissa_free(temp);
      temp=NULL;
    }

  set_borders_invalid(out);
}
