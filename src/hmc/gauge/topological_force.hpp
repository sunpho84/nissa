#ifndef _TOPOLOGICAL_FORCE_HPP
#define _TOPOLOGICAL_FORCE_HPP

#include "new_types/su3.hpp"
#include "topological_action.hpp"

namespace nissa
{
  void compute_topological_force_lx_conf(quad_su3 *F,quad_su3 *conf,topotential_pars_t *pars);
}

#endif
