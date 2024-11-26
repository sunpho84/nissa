#ifndef _CGM_INVERT_STD2EE_M2_PORTABLE_HPP
#define _CGM_INVERT_STD2EE_M2_PORTABLE_HPP

#include <dirac_operators/stD/dirac_operator_stD.hpp>
#include <inverters/templates/cgm_invert_template_threaded.hpp>
#include <new_types/rat_approx.hpp>

namespace nissa
{
  inline void inv_stD2ee_m2_cgm_portable_run_hm_up_to_comm_prec(std::vector<EvnField<color>>& chi_e,
								const EoField<quad_su3>& eo_conf,
								const std::vector<double>& poles,
								const int& niter_max,
								const double& residue,
								const EvnField<color>& pf)
  {
    cgm_invert(chi_e,
	       poles,
	       [temp=OddField<color>("temp",WITH_HALO),
		&eo_conf]
	       (EvnField<color>& out,
		const double& mass2,
		const EvnField<color>& in) mutable
	       {
		 apply_stD2ee_m2(out,eo_conf,temp,mass2,in);
	       },
	       niter_max,
	       std::vector<double>(poles.size(),residue),
	       pf);
  }
  
  inline void summ_src_and_all_inv_stD2ee_m2_cgm_portable(EvnField<color>& chi_e,
							  const EoField<quad_su3>& eo_conf,
							  const rat_approx_t& appr,
							  const int& niter_max,
							  const double& req_res,
							  const EvnField<color>& source)
  {
    std::vector<EvnField<color>> temp(appr.degree(),{"temp",WITH_HALO});
    
    inv_stD2ee_m2_cgm_portable_run_hm_up_to_comm_prec(temp,eo_conf,appr.poles,niter_max,req_res,source);
    
    chi_e=source;
    chi_e*=appr.cons;
    
    for(int iterm=0;iterm<appr.degree();iterm++)
      FOR_EACH_SITE_DEG_OF_FIELD(chi_e,
				 CAPTURE(appr,
					 TO_WRITE(chi_e),
					 w=appr.weights[iterm],
					 t=temp[iterm].getReadable()),
				 site,iDeg,
				 {
				   chi_e(site,iDeg)+=
				     w*t(site,iDeg);
				 });
  }
}

#endif
