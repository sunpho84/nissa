#ifndef _DIRAC_OPERATOR_TMCLOVQ2_128_HPP
#define _DIRAC_OPERATOR_TMCLOVQ2_128_HPP

#include "new_types/float_128.hpp"

namespace nissa
{
  void apply_tmclovQ2_128(spincolor_128 *out,quad_su3 *conf,double kappa,clover_term_t *Cl,spincolor_128 *temp,double mu,spincolor_128 *in);
}

#endif
