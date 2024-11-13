#ifndef _CG_INVERT_EVN_STD_HPP
#define _CG_INVERT_EVN_STD_HPP

#include <optional>

#include "base/field.hpp"

namespace nissa
{
  void inv_evn_stD_cg(EvnField<color>& sol,
		      const std::optional<EvnField<color>>& guess,
		      const EoField<quad_su3>& conf,
		      const double& m,
		      const int& niter,
		      const double& residue,
		      const EoField<color>& source);
}

#endif
