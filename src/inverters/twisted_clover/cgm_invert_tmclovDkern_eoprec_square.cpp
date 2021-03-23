#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "cgm_invert_tmclvDkern_eoprec_square_portable.hpp"

#include "new_types/rat_approx.hpp"

namespace nissa
{
  void summ_src_and_all_inv_stD2ee_m2_cgm(color *chi_e,quad_su3 **eo_conf,rat_approx_t *appr,int niter_max,double req_res,color *source)
  {
    summ_src_and_all_inv_stD2ee_m2_cgm_portable(chi_e,eo_conf,appr,niter_max,req_res,source);
  }
  
  void inv_stD2ee_m2_cgm_run_hm_up_to_comm_prec(color **chi_e,quad_su3 **eo_conf,double *poles,int nterms,int niter_max,double residue,color *pf)
  {
    inv_stD2ee_m2_cgm_portable_run_hm_up_to_comm_prec(chi_e,eo_conf,poles,nterms,niter_max,residue,pf);
  }
}
