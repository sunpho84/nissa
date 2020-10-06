#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/communicate.hpp"

#include "../WclovQ/dirac_operator_WclovQ.hpp"

namespace nissa
{
  void apply_WclovQ2(spincolor *out,quad_su3 *conf,double kappa,clover_term_t *Cl,spincolor *temp,spincolor *in)
  {
    apply_WclovQ(temp,conf,kappa,Cl,in);
    apply_WclovQ(out,conf,kappa,Cl,temp);
  } 
}
