#ifndef _CG_INVERT_TMCLOVD_EOPREC_HPP
#define _CG_INVERT_TMCLOVD_EOPREC_HPP

#include <optional>

#include <base/field.hpp>

namespace nissa
{
  void inv_tmclovD_cg_eoprec(LxField<spincolor>& solution_lx,
			     std::optional<OddField<spincolor>> guess_Koo,
			     const LxField<quad_su3>& conf_lx,
			     const double& kappa,
			     const LxField<clover_term_t>& Cl_lx,
			     const LxField<inv_clover_term_t>* ext_invCl_lx,
			     const double& cSW,
			     const double& mass,
			     const int& nitermax,
			     const double& targResidue,
			     const LxField<spincolor>& source_lx);
}

#endif
