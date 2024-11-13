#ifndef _DIRAC_OPERATOR_TMCLOVQ_HPP
#define _DIRAC_OPERATOR_TMCLOVQ_HPP

#include "base/field.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void apply_tmclovQ(LxField<spincolor>& out,
		     const LxField<quad_su3>& conf,
		     const double& kappa,
		     const LxField<clover_term_t>& Cl,
		     const double& mu,
		     const LxField<spincolor>& in);
}

#endif
