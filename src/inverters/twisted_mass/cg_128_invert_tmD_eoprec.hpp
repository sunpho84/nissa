#ifndef _CG_128_INVERT_TMDKERN_EOPREC_SQUARE_EOS_HPP
#define _CG_128_INVERT_TMDKERN_EOPREC_SQUARE_EOS_HPP

#include <optional>

#include "base/field.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmDkern_eoprec_square_eos_cg_128(OddField<spincolor>& sol,
					    std::optional<OddField<spincolor>> guess,
					    const EoField<quad_su3>& conf,
					    const double& kappa,
					    const double& mu,
					    const int& niter,
					    const double external_solver_residue,
					    const OddField<spincolor>& external_source);
}

#endif
