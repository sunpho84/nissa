#ifndef _CGM_32_INVERT_STD2EE_M2_HPP
#define _CGM_32_INVERT_STD2EE_M2_HPP

namespace nissa
{
  void summ_src_and_all_inv_stD2ee_m2_cgm_32(color *chi_e,quad_su3 **eo_conf,rat_approx_t *appr,int niter_max,double req_res,color *source);
  void summ_src_and_all_inv_stD2ee_m2_cgm_32(color *chi_e,quad_su3 **eo_conf,rat_approx_t *appr,int niter_max,double req_res,color *source);
  void inv_stD2ee_m2_cgm_32_run_hm_up_to_comm_prec(color **chi_e,quad_su3 **eo_conf,double *poles,int nterms,int niter_max,double residue,color *pf);
}

#endif
