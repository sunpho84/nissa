#ifndef _CG_EOPREC_TWISTED_MASS_FREE_OPERATOR_HPP
#define _CG_EOPREC_TWISTED_MASS_FREE_OPERATOR_HPP

#include <optional>

#include "base/field.hpp"
#include "free_theory_types.hpp"

namespace nissa
{
  void inv_tmDkern_eoprec_square_eos(OddField<spin0>& sol,
				     const std::optional<OddField<spin0>>& guess,
				     const tm_quark_info& qu,
				     const int& nMaxIter,
				     const double& residue,
				     const OddField<spin0>& source);

  void inv_tmD_cg_eoprec_eos(LxField<spin0>& solution_lx,
			     std::optional<OddField<spin0>> guess_Koo,
			     const tm_quark_info& qu,
			     const int& nitermax,
			     const double& residue,
			     const LxField<spin0>& source_lx);
}

#endif
