#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "cg_invert_stD2ee_m2_portable.hpp"

namespace nissa
{
  void inv_stD2ee_m2_cg(EvnField<color>& sol,
			const std::optional<EvnField<color>>& guess,
			const EoField<quad_su3>& eo_conf,
			const double& m2,
			const int& niter,
			const double& residue,
			const EvnField<color>& source)
  {
    inv_stD2ee_m2_cg_portable(sol,guess,eo_conf,m2,niter,residue,source);
  }
}
