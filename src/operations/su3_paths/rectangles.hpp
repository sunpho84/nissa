#ifndef _RECTANGLES_HPP
#define _RECTANGLES_HPP

#include "base/old_field.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void global_plaquette_and_rectangles_lx_conf(double* glb_shapes,
					       const LxField<quad_su3>& conf);
  void point_plaquette_and_rectangles_lx_conf(LxField<complex>& point_shapes,
					      const LxField<quad_su3>& conf);
  void global_plaquette_and_rectangles_lx_conf_per_timeslice(double* glb_shapes,
							     const LxField<quad_su3>& conf);
  void global_plaquette_and_rectangles_eo_conf(double *paths,eo_ptr<quad_su3> conf);
}

#endif
