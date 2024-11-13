#ifndef _CG_INVERT_STD2_M2_H
#define _CG_INVERT_STD2_M2_H

#include <optional>

#include "base/field.hpp"
#include "geometry/geometry_eo.hpp"

namespace nissa
{
  void inv_stD_cg(EoField<color>& sol,
		  const std::optional<EvnField<color>>& guess,
		  const EoField<quad_su3>& conf,
		  const double& m,
		  const int& niter,
		  const double& residue,
		  const EoField<color>& source);
  
  inline void inv_stD_cg(EoField<color>& sol,
			 const EoField<quad_su3>& conf,
			 const double& m,
			 const int& niter,
			 const double& residue,
			 const EoField<color>& source)
  {
    inv_stD_cg(sol,nullptr,conf,m,niter,residue,source);
  }
}

#endif
