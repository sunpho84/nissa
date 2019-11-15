#ifndef _PLAQUETTE_HPP
#define _PLAQUETTE_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  double global_plaquette_lx_conf(quad_su3 *conf);
  double global_plaquette_eo_conf(quad_su3 **conf);
  void global_plaquette_eo_conf(double *totplaq,quad_su3 **conf);
  void global_plaquette_lx_conf(double *totplaq,quad_su3 *conf);
}

#endif
