#ifndef _DIRAC_OPERATOR_HOVERLAP_HPP
#define _DIRAC_OPERATOR_HOVERLAP_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_H(spincolor *out,quad_su3 *conf,double kappa, spincolor *in);
}

#endif
