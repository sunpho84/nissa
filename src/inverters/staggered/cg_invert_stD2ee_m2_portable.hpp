#ifndef _CG_INVERT_STD2EE_M2_PORTABLE_HPP
#define _CG_INVERT_STD2EE_M2_PORTABLE_HPP

#include <optional>

#include <base/old_field.hpp>
#include <new_types/su3.hpp>

namespace nissa
{
  void inv_stD2ee_m2_cg_portable(EvnField<color>& sol,
				 const std::optional<EvnField<color>>& guess,
				 const EoField<quad_su3> conf,
				 const double m2,
				 const int& niter,
				 const double& residue,
				 const EvnField<color>& source);
}

#endif
