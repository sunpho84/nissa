#ifndef _MULTIPSEUDO_RHMC_STEP_HPP
#define _MULTIPSEUDO_RHMC_STEP_HPP

namespace nissa
{
  double multipseudo_rhmc_step(quad_su3 **out_conf,quad_su3 **in_conf,theory_pars_t &theory_pars,
				 hmc_evol_pars_t &simul_pars,int itraj);
}

#endif
