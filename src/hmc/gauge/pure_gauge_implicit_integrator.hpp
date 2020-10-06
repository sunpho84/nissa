#ifndef _PURE_GAUGE_IMPLICIT_INTEGRATOR_HPP
#define _PURE_GAUGE_IMPLICIT_INTEGRATOR_HPP

#include "hmc/gauge/pure_gauge_Omelyan_integrator.hpp"
#include "hmc/theory_pars.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void implicit_pure_gauge_evolver(quad_su3 *H,quad_su3 *conf,su3 **pi,su3 **phi,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *simul);
  void implicit_pure_gauge_leapfrog_evolver(quad_su3 *H,quad_su3 *conf,su3 **pi,su3 **phi,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *simul);
}

#endif
