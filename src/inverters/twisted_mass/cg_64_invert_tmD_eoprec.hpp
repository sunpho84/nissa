#ifndef _CG_64_INVERT_TMDEOIMPR_HPP
#define _CG_64_INVERT_TMDEOIMPR_HPP

#include <optional>

#include <base/field.hpp>

namespace nissa
{
  void inv_tmDkern_eoprec_square_eos_cg_64(OddField<spincolor>& sol,
					   std::optional<OddField<spincolor>> guess,
					   const EoField<quad_su3>& conf,
					   const double& kappa,
					   const double& mu,
					   const int& niter,
					   const double residue,
					   const OddField<spincolor>& source);
}

#endif
