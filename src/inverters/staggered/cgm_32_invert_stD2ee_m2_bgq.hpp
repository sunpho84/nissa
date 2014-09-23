#ifndef _CGM_SINGLE_INVERT_STD2EE_M2_BGQ_H
#define _CGM_SINGLE_INVERT_STD2EE_M2_BGQ_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void summ_src_and_all_inv_stD2ee_m2_cgm_32_bgq(bi_single_color *chi_e,bi_single_oct_su3 **eo_conf,rat_approx_t *appr,int niter_max,double req_res,bi_single_color *source);
  void inv_stD2ee_m2_cgm_32_bgq(bi_single_color** sol,bi_single_oct_su3** conf,double* shift,int nshift,int niter_max,double* ext_req_res,bi_single_color* source);
  void inv_stD2ee_m2_cgm_32_bgq_run_hm_up_to_comm_prec(bi_single_color **chi_e,bi_single_oct_su3 **eo_conf,double *poles,int nterms,int niter_max,double residue,bi_single_color *pf);  
}

#endif
