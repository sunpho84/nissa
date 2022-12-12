#ifndef _CG_INVERT_TMCLOVQ2_64_PORTABLE_HPP
#define _CG_INVERT_TMCLOVQ2_64_PORTABLE_HPP

#include <math.h>
#include <cmath>
#include <optional>

#include <base/field.hpp>
#include <dirac_operators/tmclovQ2/dirac_operator_tmclovQ2.hpp>
#include <inverters/templates/cg_invert_template_threaded.hpp>

namespace nissa
{
  inline
  void inv_tmclovQ2_cg_64_portable(LxField<spincolor>& sol,
				   std::optional<LxField<spincolor>> guess,
				   const LxField<quad_su3>& conf,
				   const double& kappa,
				   const LxField<clover_term_t>& Cl,
				   const double& mu,
				   const int& niter,
				   const double& residue,
				   const LxField<spincolor>& source)
  {
    cg_invert(sol,
	      guess,
	      [temp=LxField<spincolor>("temp",WITH_HALO),
	       &conf,
	       &kappa,
	       &Cl,
	       &mu](LxField<spincolor>& out,
		    const LxField<spincolor>& in) mutable
	      {
		apply_tmclovQ2(out,conf,kappa,Cl,temp,mu,in);
	      },
	      niter,
	      residue,
	      source);
  }
}

#endif
