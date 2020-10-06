#ifndef _RECTANGLE_STAPLE_HPP
#define _RECTANGLE_STAPLE_HPP

#include "squared_staples.hpp"

namespace nissa
{
  typedef squared_staples_t rectangular_staples_t;

  void compute_rectangular_staples_lx_conf(rectangular_staples_t *out,quad_su3 *conf,squared_staples_t *sq_staples);
  void compute_summed_rectangular_staples_lx_conf(quad_su3 *out,quad_su3 *conf,squared_staples_t *sq_staples);
}

#endif
