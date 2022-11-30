#ifndef _CG_INVERT_TMD_EOPREC_HPP
#define _CG_INVERT_TMD_EOPREC_HPP

#include <optional>

#include "base/field.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmD_cg_eoprec_prepare_source(spincolor *varphi,eo_ptr<quad_su3> conf_eos,spincolor *eq8a,spincolor *source_odd);

  void inv_tmD_cg_eoprec(LxField<spincolor>& solution_lx,
			 std::optional<OddField<spincolor>> guess_Koo,
			 const LxField<quad_su3>& conf_lx,
			 const double& kappa,
			 const double& mass,
			 const int& nitermax,
			 const double& residue,
			 const LxField<spincolor>& source_lx);
  
  void inv_tmDkern_eoprec_square_eos(spincolor *sol,spincolor *guess,eo_ptr<quad_su3> conf,double kappa,double mu,int nitermax,double residue,spincolor *source);
  void inv_tmD_cg_eoprec_almost_reco_sol(spincolor *varphi,eo_ptr<quad_su3> conf_eos,spincolor *sol_odd,spincolor *source_evn);
}

#endif
