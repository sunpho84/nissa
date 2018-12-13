#ifndef _CGM_INVERT_OVERLAP_KERNEL2_HPP
#define _CGM_INVERT_OVERLAP_KERNEL2_HPP

#include "new_types/rat_approx.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void inv_overlap_kernel2_cgm(spincolor **sol,quad_su3 *conf,double M,double *shift,int nshift,int niter_max,double *ext_req_res,spincolor *source);
  void summ_src_and_all_inv_overlap_kernel2_cgm(spincolor *sol,quad_su3 *conf,double M,rat_approx_t *appr,int niter_max,double req_res,spincolor *source);
}

#endif
