#ifndef _DIRAC_OPERATOR_TMCLOVQ2_128_HPP
#define _DIRAC_OPERATOR_TMCLOVQ2_128_HPP

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void apply_tmclovQ2_128(spincolor_128 *out,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,spincolor_128 *temp,double mu,spincolor_128 *in);
  void apply_tmclovQ2_128(spincolor_128 *out,quad_su3 *conf,double kappa,opt_as2t_su3 *Cl,spincolor_128 *temp,double mu,spincolor_128 *in);
}

#endif
