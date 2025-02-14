#ifndef _CG_INVERT_STD2EE_M2_HPP
#define _CG_INVERT_STD2EE_M2_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <optional>

#include "base/field.hpp"
#include "inverters/templates/cg_invert_template_threaded.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  template <typename F>
  void inv_stD2ee_m2_cg(EvnField<Color<F>>& sol,
			const std::optional<EvnField<Color<F>>>& guess,
			const EoField<QuadSu3<F>>& conf,
			const double& m2,
			const int& niter,
			const double& residue,
			const EvnField<Color<F>>& source)
  {
    std::function<void(EvnField<Color<F>>&,
		       const EvnField<Color<F>>&)> f=
      [temp=OddField<Color<F>>("temp",WITH_HALO),
       &conf,
       &m2]
      (EvnField<Color<F>>& out,
       const EvnField<Color<F>>& in) mutable
      {
	apply_stD2ee_m2(out,conf,temp,m2,in);
      };
    
    cg_invert(sol,
	      guess,
	      f,
	      niter,
	      residue,
	      source);
  }
}

#endif
