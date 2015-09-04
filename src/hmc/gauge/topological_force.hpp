#ifndef _TOPOLOGICAL_FORCE_H
#define _TOPOLOGICAL_FORCE_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void compute_topological_force_lx_conf(quad_su3 *F,quad_su3 *conf,topotential_pars_t *pars,bool phase_pres);
}

#endif
