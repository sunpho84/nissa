#ifndef _CGM_INVERT_STD2EE_M2_PORTABLE_HPP
#define _CGM_INVERT_STD2EE_M2_PORTABLE_HPP

#include "geometry/geometry_eo.hpp"
#include "new_types/rat_approx.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void summ_src_and_all_inv_stD2ee_m2_cgm_portable(color *chi_e,eo_ptr<quad_su3> eo_conf,rat_approx_t *appr,int niter_max,double req_res,color *source);
  void inv_stD2ee_m2_cgm_portable_run_hm_up_to_comm_prec(eo_ptr<color> chi_e,eo_ptr<quad_su3> eo_conf,double *poles,int nterms,int niter_max,double residue,color *pf);
}

#endif
