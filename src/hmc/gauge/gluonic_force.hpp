#ifndef _GLUONIC_FORCE_HPP
#define _GLUONIC_FORCE_HPP

#include "hmc/theory_pars.hpp"

namespace nissa
{
  void gluonic_force_finish_computation(quad_su3*F,quad_su3*conf);
  void compute_gluonic_force_lx_conf(quad_su3*F,quad_su3*conf,theory_pars_t*physics);
}

#endif
