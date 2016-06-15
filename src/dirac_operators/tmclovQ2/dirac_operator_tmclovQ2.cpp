#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "dirac_operators/tmclovQ/dirac_operator_tmclovQ.hpp"
#include "new_types/su3.hpp"

//Apply the Q+Q- operator to a spincolor

namespace nissa
{
  void apply_tmclovQ2(spincolor *out,quad_su3 *conf,double kappa,clover_term_t *Cl,spincolor *temp,double mu,spincolor *in)
  {
    apply_tmclovQ(temp,conf,kappa,Cl,+mu,in);
    apply_tmclovQ(out,conf,kappa,Cl,-mu,temp);
  }
  
  void apply_tmclovQ2_m2(spincolor *out,quad_su3 *conf,double kappa,clover_term_t *Cl,spincolor *temp,double mu,spincolor *in)
  {apply_tmclovQ2(out,conf,kappa,Cl,temp,sqrt(mu),in);}
}
