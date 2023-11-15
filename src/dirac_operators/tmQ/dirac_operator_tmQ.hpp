#ifndef _DIRAC_OPERATOR_TMQ_HPP
#define _DIRAC_OPERATOR_TMQ_HPP

#include "base/old_field.hpp"

namespace nissa
{
  void apply_tmQ(LxField<spincolor>& out,
		 const LxField<quad_su3>& conf,
		 const double& kappa,
		 const double& mu,
		 const LxField<spincolor>& in);
    
  void apply_tmQ_v1(spincolor *out,quad_su3 *conf,double kappa,double mu,spincolor *in);
}

#endif
