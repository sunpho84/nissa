#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "dirac_operators/tmclovQ/dirac_operator_tmclovQ_128.hpp"

//Apply the Q+Q- operator to a spincolor

namespace nissa
{
  void apply_tmclovQ2_128(spincolor_128 *out,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,spincolor_128 *temp,double mu,spincolor_128 *in)
  {
    apply_tmclovQ_128(temp,conf,kappa,csw,Pmunu,+mu,in);
    apply_tmclovQ_128(out,conf,kappa,csw,Pmunu,-mu,temp);
  }

  void apply_tmclovQ2_128(spincolor_128 *out,quad_su3 *conf,double kappa,opt_as2t_su3 *Cl,spincolor_128 *temp,double mu,spincolor_128 *in)
  {
    apply_tmclovQ_128(temp,conf,kappa,Cl,+mu,in);
    apply_tmclovQ_128(out,conf,kappa,Cl,-mu,temp);
  }
}
