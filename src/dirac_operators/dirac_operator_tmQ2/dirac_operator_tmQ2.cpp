#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <math.h>

#include "../dirac_operator_tmQ/dirac_operator_tmQ.h"
#include "../dirac_operator_tmQ_left/dirac_operator_tmQ_left.h"

#include "../../communicate/communicate.h"
#include "../../base/global_variables.h"
#include "../../base/thread_macros.h"
#include "../../base/vectors.h"
#include "../../new_types/new_types_definitions.h"
#ifdef USE_THREADS
 #include "../../routines/thread.h"
#endif

//Apply the Q+Q- operator to a spincolor
THREADABLE_FUNCTION_7ARG(apply_tmQ2_RL, spincolor*,out, quad_su3*,conf, double,kappa, spincolor*,ext_temp, int,RL, double,mu, spincolor*,in)
{
  spincolor *temp=ext_temp;
  if(temp==NULL) temp=nissa_malloc("tempQ",loc_vol+bord_vol,spincolor);
  
  if(RL==0) apply_tmQ(temp,conf,kappa,+mu,in);
  else apply_tmQ_left(temp,conf,kappa,+mu,in);
  
  if(RL==0) apply_tmQ(out,conf,kappa,-mu,temp);
  else apply_tmQ_left(out,conf,kappa,-mu,temp);

  if(ext_temp==NULL) nissa_free(temp);
}}

//wrappers
void apply_tmQ2_m2_RL(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,int RL,double m2,spincolor *in)
{apply_tmQ2_RL(out,conf,kappa,temp,RL,sqrt(m2),in);}

void apply_tmQ2(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,double mu,spincolor *in)
{apply_tmQ2_RL(out,conf,kappa,temp,0,mu,in);}

void apply_tmQ2_left(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,double mu,spincolor *in)
{apply_tmQ2_RL(out,conf,kappa,temp,1,mu,in);}
