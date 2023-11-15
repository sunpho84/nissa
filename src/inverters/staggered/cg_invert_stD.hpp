#ifndef _CG_INVERT_STD2_M2_H
#define _CG_INVERT_STD2_M2_H

#include <optional>

#include "base/old_field.hpp"
#include "geometry/geometry_eo.hpp"

namespace nissa
{
  void inv_stD_cg(OldEoField<color0>& sol,
		  const std::optional<EvnField<color0>>& guess,
		  const OldEoField<quad_su3>& conf,
		  const double& m,
		  const int& niter,
		  const double& residue,
		  const OldEoField<color0>& source);
  
  inline void inv_stD_cg(OldEoField<color0>& sol,
			 const OldEoField<quad_su3>& conf,
			 const double& m,
			 const int& niter,
			 const double& residue,
			 const OldEoField<color0>& source)
  {
    inv_stD_cg(sol,nullptr,conf,m,niter,residue,source);
  }
}

#endif
