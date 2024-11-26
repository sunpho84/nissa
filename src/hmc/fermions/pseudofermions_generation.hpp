#ifndef _PSEUDOFERMIONS_GENERATION_HPP
#define _PSEUDOFERMIONS_GENERATION_HPP

#include "hmc/multipseudo/multipseudo_rhmc_step.hpp"

namespace nissa
{
  double generate_pseudofermions(std::vector<std::vector<pseudofermion_t>>& pf,
				 EoField<quad_su3>& conf,
				 const theory_pars_t& theory_pars,
				 const hmc_evol_pars_t& simul_pars,
				 const std::vector<rat_approx_t>& rat_appr);
}

#endif
