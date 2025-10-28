#ifndef _DIRAC_OPERATOR_TMQ2_128_HPP
#define _DIRAC_OPERATOR_TMQ2_128_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <new_types/float_128.hpp>

namespace nissa
{
  void apply_tmQ2_128(SpinColor128 *out,quad_su3 *conf,double kappa,SpinColor128 *temp,double mu,SpinColor128 *in);
}

#endif
