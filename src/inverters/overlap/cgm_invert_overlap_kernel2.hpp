#ifndef _CGM_INVERT_OVERLAP_KERNEL2_HPP
#define _CGM_INVERT_OVERLAP_KERNEL2_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_overlap_kernel2_cgm(spincolor **sol,quad_su3 *conf,double M, int niter_max,double *req_res,spincolor *source);
}

#endif
