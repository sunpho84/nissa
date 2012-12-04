#include "../../new_types/new_types_definitions.h"
#include "../../base/vectors.h"
#include "../../base/communicate.h"
#include "../../base/global_variables.h"

#include "../dirac_operator_WclovQ/dirac_operator_WclovQ.h"

void apply_WclovQ2(spincolor *out,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,spincolor *ext_temp,spincolor *in)
{
  spincolor *temp=(ext_temp!=NULL)?ext_temp:nissa_malloc("temp",loc_vol+bord_vol,spincolor);
  
  apply_WclovQ(temp,conf,kappa,csw,Pmunu,in);
  apply_WclovQ(out,conf,kappa,csw,Pmunu,temp);
  
  if(ext_temp==NULL) nissa_free(temp);
}
