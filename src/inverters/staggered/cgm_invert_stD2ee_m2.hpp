#ifndef _CGM_INVERT_STD2EE_M2_HPP
#define _CGM_INVERT_STD2EE_M2_HPP

#include "cgm_invert_stD2ee_m2_portable.hpp"

namespace nissa
{
  inline
  void summ_src_and_all_inv_stD2ee_m2_cgm(EvnField<color>& chi_e,
					  const EoField<quad_su3>& eo_conf,
					  const rat_approx_t& appr,
					  const int& niter_max,
					  const double& req_res,
					  const EvnField<color>& source)
  {
    summ_src_and_all_inv_stD2ee_m2_cgm_portable(chi_e,eo_conf,appr,niter_max,req_res,source);
  }
  
  inline
  void inv_stD2ee_m2_cgm_run_hm_up_to_comm_prec(std::vector<EvnField<color>>& chi_e,
						const EoField<quad_su3>& eo_conf,
						const std::vector<double>& poles,
						const int& niter_max,
						const double& residue,
						const EvnField<color>& pf)
  {
    inv_stD2ee_m2_cgm_portable_run_hm_up_to_comm_prec(chi_e,eo_conf,poles,niter_max,residue,pf);
  }
}

#endif
