#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/new_types_definitions.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "base/global_variables.hpp"

#include "dirac_operator_WclovQ/dirac_operator_WclovQ.hpp"

namespace nissa
{
  void apply_WclovQ2(spincolor *out,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,spincolor *temp,spincolor *in)
  {
    apply_WclovQ(temp,conf,kappa,csw,Pmunu,in);
    apply_WclovQ(out,conf,kappa,csw,Pmunu,temp);
  } 
}
