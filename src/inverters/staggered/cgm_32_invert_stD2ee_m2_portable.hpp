#ifndef _CGM_SINGLE_INVERT_STD2EE_M2_PORTABLE_HPP
#define _CGM_SINGLE_INVERT_STD2EE_M2_PORTABLE_HPP

#include "new_types/rat_approx.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void summ_src_and_all_inv_stD2ee_m2_cgm_32_portable(single_color *chi_e,single_quad_su3 **eo_conf,rat_approx_t *appr,int niter_max,double req_res,single_color *source);
  void inv_stD2ee_m2_cgm_32_portable(single_color** sol,single_quad_su3** conf,double* shift,int nshift,int niter_max,double* ext_req_res,single_color* source);
  void inv_stD2ee_m2_cgm_32_portable_run_hm_up_to_comm_prec(single_color **chi_e,single_quad_su3 **eo_conf,double *poles,int nterms,int niter_max,double residue,single_color *pf);
}

#endif
