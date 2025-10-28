#ifndef _CG_INVERT_TMCLOVQ2_HPP
#define _CG_INVERT_TMCLOVQ2_HPP

#include "cg_64_invert_tmclovQ2.hpp"

namespace nissa
{
  //switch 64 and 128
  inline void inv_tmclovQ2_cg(LxField<spincolor>& sol,
			      std::optional<LxField<spincolor>> guess,
			      const LxField<quad_su3>& conf,
			      const double& kappa,
			      const LxField<clover_term_t>& Cl,
			      const double& mu,
			      const int& niter,
			      const double& residue,
			      const LxField<spincolor>& source)
  {
    if(use_128_bit_precision)
      CRASH("reimplement");//inv_tmclovQ2_cg_128(sol,guess,conf,kappa,Cl,m,niter,residue,source);
    else inv_tmclovQ2_cg_64(sol,guess,conf,kappa,Cl,mu,niter,residue,source);
  }
}

#endif
