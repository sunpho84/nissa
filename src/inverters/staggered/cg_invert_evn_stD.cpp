#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <optional>

#include "base/field.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "inverters/staggered/cg_invert_stD2ee_m2.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void inv_evn_stD_cg(EvnField<color>& sol,
		      const std::optional<EvnField<color>>& guess,
		      const EoField<quad_su3>& conf,
		      const double& m,
		      const int& niter,
		      const double& residue,
		      const EoField<color>& source)
  {
    //apply the dagger ...
    EvnField<color> temp("temp",WITH_HALO);
    evn_apply_stD_dag(temp,conf,m,source);
    
    //and invert the DD^+
    inv_stD2ee_m2_cg(sol,guess,conf,m*m,niter,residue,temp);
  }
}
