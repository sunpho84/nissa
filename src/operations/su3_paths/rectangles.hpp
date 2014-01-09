#ifndef _RECTANGLES_H
#define _RECTANGLES_H

namespace nissa
{
  void point_plaquette_and_rectangles_lx_conf(double *paths,quad_su3 *conf);
  void global_plaquette_and_rectangles_lx_conf(double *paths,quad_su3 *conf);
  void global_plaquette_and_rectangles_lx_conf_per_timeslice(double *paths,quad_su3 *conf);
  void global_plaquette_and_rectangles_eo_conf(double *paths,quad_su3 **conf);
}

#endif
