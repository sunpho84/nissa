#ifndef _OMELYAN_INTEGRATOR_HPP
#define _OMELYAN_INTEGRATOR_HPP

#include "multipseudo_rhmc_step.hpp"

namespace nissa
{
  void Omelyan_integrator(EoField<quad_su3>& H,
			  EoField<quad_su3>& conf,
			  std::vector<std::vector<pseudofermion_t>>& pf,
			  theory_pars_t& theory_pars,
			  hmc_evol_pars_t& simul_pars,
			  std::vector<rat_approx_t>& rat_appr);
}

#endif
