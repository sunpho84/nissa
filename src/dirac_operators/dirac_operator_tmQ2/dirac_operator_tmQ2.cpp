#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <math.h>

#include "../../new_types/new_types_definitions.h"
#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/vectors.h"
#include "../dirac_operator_tmQ/dirac_operator_tmQ.h"
#include "../dirac_operator_tmQ_left/dirac_operator_tmQ_left.h"

//WORKAROUND
int all=0;
  
//Apply the Q+Q- operator to a spincolor
void apply_tmQ2_RL(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,int RL,double mu,spincolor *in)
{
#pragma omp single
  if(temp==NULL)
    {
      temp=nissa_malloc("tempQ",loc_vol+bord_vol,spincolor);
      all=1;
    }
  
  if(RL==0) apply_tmQ(temp,conf,kappa,+mu,in);
  else apply_tmQ_left(temp,conf,kappa,+mu,in);
  
  if(RL==0) apply_tmQ(out,conf,kappa,-mu,temp);
  else apply_tmQ_left(out,conf,kappa,-mu,temp);

#pragma omp single
  if(all==1)
    {
      nissa_free(temp);
      temp=NULL;
    }
}

//wrappers
void apply_tmQ2_m2_RL(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,int RL,double m2,spincolor *in)
{apply_tmQ2_RL(out,conf,kappa,temp,RL,sqrt(m2),in);}

void apply_tmQ2(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,double mu,spincolor *in)
{apply_tmQ2_RL(out,conf,kappa,temp,0,mu,in);}

void apply_tmQ2_left(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,double mu,spincolor *in)
{apply_tmQ2_RL(out,conf,kappa,temp,1,mu,in);}
