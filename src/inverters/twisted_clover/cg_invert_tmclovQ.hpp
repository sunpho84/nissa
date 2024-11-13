#ifndef _CG_INVERT_TMCLOVQ_HPP
#define _CG_INVERT_TMCLOVQ_HPP

#include "cg_invert_tmclovQ2.hpp"

namespace nissa
{
  inline void inv_tmclovQ_cg(LxField<spincolor>& sol,
			     std::optional<LxField<spincolor>>& guess,
			     const LxField<quad_su3>& conf,
			     const double& kappa,
			     const LxField<clover_term_t>& Cl,
			     const double& mu,
			     const int& niter,
			     const double& residue,
			     const LxField<spincolor>& source)
  {
    inv_tmclovQ2_cg(sol,guess,conf,kappa,Cl,mu,niter,residue,source);
    LxField<spincolor> temp("temp",WITH_HALO);
    
    //remove the "wrong r"
    temp=sol;
    apply_tmclovQ(sol,conf,kappa,Cl,-mu,temp);
  }
}

#endif
