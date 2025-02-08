#ifndef _DIRAC_OPERATOR_TMQ_128_H
#define _DIRAC_OPERATOR_TMQ_128_H

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include "new_types/float_128.hpp"

namespace nissa
{
  void apply_tmQ_128(SpinColor128 *out,quad_su3 *conf,double kappa,double mu,SpinColor128 *in);
}

#endif
