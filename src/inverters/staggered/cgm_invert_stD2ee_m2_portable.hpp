#ifndef _CGM_INVERT_STD2EE_M2_PORTABLE_H
#define _CGM_INVERT_STD2EE_M2_PORTABLE_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void summ_src_and_all_inv_stD2ee_m2_cgm_portable(color *chi_e,quad_su3 **eo_conf,rat_approx_t *appr,int niter_max,double req_res,color *source);
  void inv_stD2ee_m2_cgm_portable_run_hm_up_to_comm_prec(color **chi_e,quad_su3 **eo_conf,double *poles,int nterms,int niter_max,double residue,color *pf);
}

#endif
