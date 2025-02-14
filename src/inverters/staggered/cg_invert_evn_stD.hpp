#ifndef _CG_INVERT_EVN_STD_HPP
#define _CG_INVERT_EVN_STD_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <optional>

#include "base/field.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "inverters/staggered/cg_invert_stD2ee_m2.hpp"

namespace nissa
{
  template <typename F>
  void inv_evn_stD_cg(EvnField<Color<F>>& sol,
		      const std::optional<EvnField<Color<F>>>& guess,
		      const EoField<QuadSu3<F>>& conf,
		      const double& m,
		      const int& niter,
		      const double& residue,
		      const EoField<Color<F>>& source)
  {
    //apply the dagger ...
    EvnField<Color<F>> temp("temp",WITH_HALO);
    evn_apply_stD_dag(temp,conf,m,source);
    
    //and invert the DD^+
    inv_stD2ee_m2_cg(sol,guess,conf,m*m,niter,residue,temp);
  }
}

#endif
