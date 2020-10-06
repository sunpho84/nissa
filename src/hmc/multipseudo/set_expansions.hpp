#ifndef _SET_EXPANSIONS_HPP
#define _SET_EXPANSIONS_HPP

#include "new_types/rat_approx.hpp"
#include "hmc/multipseudo/multipseudo_rhmc_step.hpp"

namespace nissa
{
  void set_expansions(std::vector<rat_approx_t> *rat_appr,quad_su3 **eo_conf,theory_pars_t *theory_pars,hmc_evol_pars_t *evol_pars);
}

#endif
