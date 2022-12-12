#ifndef _CG_64_INVERT_TMCLOVD_EOPREC_HPP
#define _CG_64_INVERT_TMCLOVD_EOPREC_HPP

#include <optional>

#include "base/field.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmclovDkern_eoprec_square_eos_cg_64(OddField<spincolor>& sol,
					       std::optional<OddField<spincolor>> guess,
					       const EoField<quad_su3>& conf,
					       const double& kappa,
					       const OddField<clover_term_t>& Cl_odd,
					       const EvnField<inv_clover_term_t>& invCl_evn,
					       const double& mu,
					       const int& niter,
					       const double residue,
					       const OddField<spincolor>& source);
}

#endif
