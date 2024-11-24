#ifndef _THEORY_ACTION_HPP
#define _THEORY_ACTION_HPP

#include "multipseudo_rhmc_step.hpp"
#include "new_types/rat_approx.hpp"

namespace nissa
{
  double compute_quark_action(const EoField<quad_su3>& eo_conf,
			      const std::vector<EoField<quad_u1>*>& u1b,
			      const std::vector<EvnField<pseudofermion_t>*>& pf,
			      const std::vector<quark_content_t>& quark_content,
			      const hmc_evol_pars_t& simul_pars,
			      const std::vector<rat_approx_t>& rat_appr);
  
  double full_theory_action(const EoField<quad_su3>& eo_conf,
			    const EoField<quad_su3>& sme_conf,
			    const EoField<quad_su3>& H,
			    const std::vector<std::vector<pseudofermion_t>>& pf,
			    const theory_pars_t& theory_pars,
			    const hmc_evol_pars_t& simul_pars,
			    const std::vector<rat_approx_t>& rat_appr,
			    const double external_quark_action=-1);
}

#endif
