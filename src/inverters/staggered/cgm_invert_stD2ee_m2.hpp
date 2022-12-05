#ifndef _CGM_INVERT_STD2EE_M2_HPP
#define _CGM_INVERT_STD2EE_M2_HPP

#include <base/field.hpp>
#include <new_types/su3.hpp>

namespace nissa
{
  void summ_src_and_all_inv_stD2ee_m2_cgm(color *chi_e,eo_ptr<quad_su3> eo_conf,rat_approx_t *appr,int niter_max,double req_res,color *source);
  
  void inv_stD2ee_m2_cgm_run_hm_up_to_comm_prec(std::vector<EvnField<color>>& chi_e,
						const EoField<quad_su3>& eo_conf,
						const std::vector<double>& poles,
						const int& niter_max,
						const double& residue,
						const EvnField<color>& pf);
}

#endif
