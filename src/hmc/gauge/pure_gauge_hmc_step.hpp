#ifndef _PURE_GAUGE_HMC_STEP_HPP
#define _PURE_GAUGE_HMC_STEP_HPP

#include "new_types/rat_approx.hpp"
#include "pure_gauge_Omelyan_integrator.hpp"

namespace nissa
{
  double pure_gauge_hmc_step(quad_su3 *out_conf,quad_su3 *in_conf,theory_pars_t &theory_pars,
			     pure_gauge_evol_pars_t &simul_pars,rat_approx_t *ra,int itraj);
}

#endif
