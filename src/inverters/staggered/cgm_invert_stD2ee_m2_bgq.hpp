#ifndef _CGM_INVERT_STD2EE_M2_BGQ_HPP
#define _CGM_INVERT_STD2EE_M2_BGQ_HPP

#include "new_types/rat_approx.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void summ_src_and_all_inv_stD2ee_m2_cgm_bgq(bi_color *chi_e,bi_oct_su3 **eo_conf,rat_approx_t *appr,int niter_max,double req_res,bi_color *source);
  void inv_stD2ee_m2_cgm_bgq_run_hm_up_to_comm_prec(bi_color **chi_e,bi_oct_su3 **eo_conf,double *poles,int nterms,int niter_max,double residue,bi_color *pf);
}

#endif
