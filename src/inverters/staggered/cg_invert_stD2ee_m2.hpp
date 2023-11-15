#ifndef _CG_INVERT_STD2EE_M2_HPP
#define _CG_INVERT_STD2EE_M2_HPP

#include <optional>

#include "base/old_field.hpp"

namespace nissa
{
  void inv_stD2ee_m2_cg(EvnField<color>& sol,
			const std::optional<EvnField<color>>& guess,
			const EoField<quad_su3>& eo_conf,
			const double& m2,
			const int& niter,
			const double& residue,
			const EvnField<color>& source);
}

#endif
