#ifndef _OMELYAN_INTEGRATOR_HPP
#define _OMELYAN_INTEGRATOR_HPP

namespace nissa
{
  void Omelyan_integrator(quad_su3 **H,quad_su3 **conf,pseudofermion_t *pf,theory_pars_t *theory_pars,hmc_evol_pars_t *simul);
}

#endif
