#ifndef _DIRAC_OPERATOR_TMCLOVQ_128_HPP
#define _DIRAC_OPERATOR_TMCLOVQ_128_HPP

#include "new_types/float_128.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void apply_tmclovQ_128(spincolor_128 *out,quad_su3 *conf,double kappa,clover_term_t *Cl,double mu,spincolor_128 *in);
}

#endif
