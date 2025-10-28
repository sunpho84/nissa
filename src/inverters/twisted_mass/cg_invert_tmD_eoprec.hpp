#ifndef _CG_INVERT_TMD_EOPREC_HPP
#define _CG_INVERT_TMD_EOPREC_HPP

#include <optional>

#include <base/field.hpp>
#include <dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_portable.hpp>
#include <geometry/geometry_eo.hpp>
#include <new_types/su3.hpp>

namespace nissa
{
  //Prepare the source according to Equation (8.b)
  void inv_tmD_cg_eoprec_prepare_source(OddField<spincolor>& varphi,
					const EoField<quad_su3>& conf_eos,
					const EvnField<spincolor>& eq8a,
					const OddField<spincolor>& source_odd);
  
  void inv_tmD_cg_eoprec(LxField<spincolor>& solution_lx,
			 std::optional<OddField<spincolor>> guess_Koo,
			 const LxField<quad_su3>& conf_lx,
			 const double& kappa,
			 const double& mass,
			 const int& nitermax,
			 const double& residue,
			 const LxField<spincolor>& source_lx);
  
  void inv_tmD_cg_eoprec_almost_reco_sol(EvnField<spincolor>& varphi,
					 const EoField<quad_su3>& conf_eos,
					 const OddField<spincolor>& sol_odd,
					 const EvnField<spincolor>& source_evn);
}

#endif
