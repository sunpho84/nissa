#ifndef _CGM_INVERT_MFACC_HPP
#define _CGM_INVERT_MFACC_HPP

#include "new_types/rat_approx.hpp"

namespace nissa
{
  void summ_src_and_all_inv_MFACC_cgm(su3 *sol,quad_su3 *conf,double kappa,rat_approx_t *appr,int niter_max,double req_res,su3 *source);
  void summ_src_and_all_inv_MFACC_cgm(quad_su3 *sol,quad_su3 *conf,double kappa,rat_approx_t *appr,int niter_max,double req_res,quad_su3 *source);
}

#endif
