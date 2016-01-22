#ifndef _DIRAC_OPERATOR_TMCLOVQ_128_HPP
#define _DIRAC_OPERATOR_TMCLOVQ_128_HPP

#include "new_types/float_128.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void apply_tmclovQ_128(spincolor_128 *out,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double mu,spincolor_128 *in);
  void apply_tmclovQ_128(spincolor_128 *out,quad_su3 *conf,double kappa,opt_as2t_su3 *Cl,double mu,spincolor_128 *in);
}

#endif
