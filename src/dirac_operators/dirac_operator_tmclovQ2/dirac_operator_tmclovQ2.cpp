#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <math.h>

#include "new_types/new_types_definitions.h"
#include "base/global_variables.h"
#include "communicate/communicate.h"
#include "base/vectors.h"
#include "dirac_operator_tmclovQ/dirac_operator_tmclovQ.h"

//Apply the Q+Q- operator to a spincolor

void apply_tmclovQ2(spincolor *out,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,spincolor *temp,double mu,spincolor *in)
{
  apply_tmclovQ(temp,conf,kappa,csw,Pmunu,+mu,in);
  apply_tmclovQ(out,conf,kappa,csw,Pmunu,-mu,temp);
}

void apply_tmclovQ2_m2(spincolor *out,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,spincolor *temp,double mu,spincolor *in)
{apply_tmclovQ2(out,conf,kappa,csw,Pmunu,temp,sqrt(mu),in);}
