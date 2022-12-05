#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/field.hpp>
#include <dirac_operators/stD/dirac_operator_stD.hpp>
#include <inverters/templates/cgm_invert_template_threaded.hpp>
#include <new_types/su3.hpp>

namespace nissa
{
  void inv_stD2ee_m2_cgm_portable_run_hm_up_to_comm_prec(std::vector<EvnField<color>>& chi_e,
							 const EoField<quad_su3> eo_conf,
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
}
