#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <optional>

#include <base/field.hpp>
#include <inverters/templates/cg_invert_template_threaded.hpp>
#include <dirac_operators/stD/dirac_operator_stD.hpp>

namespace nissa
{
  void inv_stD2ee_m2_cg_portable(EvnField<color>& sol,
				 const std::optional<EvnField<color>>& guess,
				 const EoField<quad_su3> conf,
				 const double m2,
				 const int& niter,
				 const double& residue,
				 const EvnField<color>& source)
  {
    std::function<void(EvnField<color>&,
      const EvnField<color>&)> f=
      [temp=OddField<color>("temp",WITH_HALO),
       &conf,
       &m2]
      (EvnField<color>& out,
		const EvnField<color>& in) mutable
      {
	apply_stD2ee_m2(out,conf,temp,m2,in);
      };
    
    cg_invert(sol,
	       guess,
	      f,
	       niter,
	       residue,
	       source);
  }
}
