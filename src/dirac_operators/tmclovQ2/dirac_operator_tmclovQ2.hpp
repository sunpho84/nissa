#ifndef _DIRAC_OPERATOR_TMCLOVQ2_HPP
#define _DIRAC_OPERATOR_TMCLOVQ2_HPP

#include "base/old_field.hpp"
#include <dirac_operators/tmclovQ/dirac_operator_tmclovQ.hpp>

namespace nissa
{
  inline void apply_tmclovQ2(LxField<spincolor>& out,
			     const LxField<quad_su3>& conf,
			     const double& kappa,
			     const LxField<clover_term_t>& Cl,
			     LxField<spincolor>& temp,
			     const double& mu,
			     const LxField<spincolor>& in)
  {
    apply_tmclovQ(temp,conf,kappa,Cl,+mu,in);
    apply_tmclovQ(out,conf,kappa,Cl,-mu,temp);
  }
  
  inline void apply_tmclovQ2_m2(LxField<spincolor>& out,
				const LxField<quad_su3>& conf,
				const double& kappa,
				const LxField<clover_term_t>& Cl,
				LxField<spincolor>& temp,
				const double& mu2,
				const LxField<spincolor>& in)
  {
    apply_tmclovQ2(out,conf,kappa,Cl,temp,sqrt(mu2),in);
  }
}

#endif
