#ifndef _CG_INVERT_TMCLOVQ2_64_HPP
#define _CG_INVERT_TMCLOVQ2_64_HPP

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"

#include "cg_64_invert_tmclovQ2_portable.hpp"

namespace nissa
{
  inline void inv_tmclovQ2_cg_64(LxField<spincolor>& sol,
				 std::optional<LxField<spincolor>> guess,
				 const LxField<quad_su3>& conf,
				 const double& kappa,
				 const LxField<clover_term_t>& Cl,
				 const double& mu,
				 const int& niter,
				 const double& residue,
				 const LxField<spincolor>& source)
  {
    inv_tmclovQ2_cg_64_portable(sol,guess,conf,kappa,Cl,mu,niter,residue,source);
  }
}

#endif
