#ifndef _DIRAC_OPERATOR_OVERLAP_KERNEL2_HPP
#define _DIRAC_OPERATOR_OVERLAP_KERNEL2_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_overlap_kernel2(spincolor* out,quad_su3* conf,double M,spincolor* ext_temp,double  diag_coeff,spincolor* in);
}

#endif
