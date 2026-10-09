#ifndef _CG_INVERT_TMCLOVD_EOPREC_HPP
#define _CG_INVERT_TMCLOVD_EOPREC_HPP

#include <base/field.hpp>
#include <dirac_operators/tmQ/dirac_operator_tmQ.hpp>

namespace nissa
{
  void inv_tmclovD_cg_eoprec(LxField<spincolor>& solution_lx,
			     std::optional<OddField<spincolor>> guess_Koo,
			     const LxField<quad_su3>& conf_lx,
			     const double& kappa,
			     const AnisDopPars& anisDopPars,
			     const double& cSW,
			     const double& mass,
			     const int& nitermax,
			     const double& targResidue,
			     const LxField<spincolor>& source_lx);
}

#endif
