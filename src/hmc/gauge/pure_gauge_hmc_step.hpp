#ifndef _PURE_GAUGE_HMC_STEP_HPP
#define _PURE_GAUGE_HMC_STEP_HPP

namespace nissa
{
  double pure_gauge_hmc_step(quad_su3 *out_conf,quad_su3 *in_conf,theory_pars_t &theory_pars,
			     pure_gauge_evol_pars_t &simul_pars,int itraj);
}

#endif
