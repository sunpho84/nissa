#ifndef _WILSON_FORCE_HPP
#define _WILSON_FORCE_HPP

#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void Wilson_force_eo_conf(eo_ptr<quad_su3> F,eo_ptr<quad_su3> eo_conf,double beta);
  void Wilson_force_lx_conf(quad_su3 *F,quad_su3 *lx_conf,double beta);
}

#endif
