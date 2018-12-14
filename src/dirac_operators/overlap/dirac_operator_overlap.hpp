#ifndef _DIRAC_OPERATOR_OVERLAP_KERNEL2_HPP
#define _DIRAC_OPERATOR_OVERLAP_KERNEL2_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_overlap(spincolor *out, quad_su3 *conf, double M, double minerr,spincolor *in);
}

#endif
