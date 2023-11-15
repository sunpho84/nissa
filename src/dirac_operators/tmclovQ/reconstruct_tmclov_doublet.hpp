#ifndef _RECONSTRUCT_TMCLOV_DOUBLET_HPP
#define _RECONSTRUCT_TMCLOV_DOUBLET_HPP

#include "base/old_field.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void reconstruct_tmclov_doublet(LxField<spincolor>& outminus,
				  LxField<spincolor>& outplus,
				  const LxField<quad_su3>& conf,
				  const double& kappa,
				  const LxField<clover_term_t>& Cl,
				  const double& mu,
				  const LxField<spincolor>& in);
}

#endif
