#ifndef _RECTANGLES_HPP
#define _RECTANGLES_HPP

#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void point_plaquette_and_rectangles_lx_conf(double *paths,quad_su3 *conf);
  void global_plaquette_and_rectangles_lx_conf(double *paths,quad_su3 *conf);
  void global_plaquette_and_rectangles_lx_conf_per_timeslice(double *paths,quad_su3 *conf);
  void global_plaquette_and_rectangles_eo_conf(double *paths,eo_ptr<quad_su3> conf);
}

#endif
