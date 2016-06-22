#ifndef _DIRAC_OPERATOR_TMQ2_128_HPP
#define _DIRAC_OPERATOR_TMQ2_128_HPP

#include "new_types/float_128.hpp"

namespace nissa
{
  void apply_tmQ2_RL_128(spincolor_128 *out,quad_su3 *conf,double kappa,spincolor_128 *temp,int RL,double mu,spincolor_128 *in);
}

#endif
