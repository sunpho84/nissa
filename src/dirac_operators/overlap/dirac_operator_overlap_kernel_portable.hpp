#ifndef _DIRAC_OPERATOR_OVERLAP_KERNEL_PORTABLE_HPP
#define _DIRAC_OPERATOR_OVERLAP_KERNEL_PORTABLE_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_overlap_kernel(spincolor *out,quad_su3 *conf,double kappa,spincolor *in);
}

#endif
