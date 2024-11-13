#ifndef _MOMENTA_ACTION_HPP
#define _MOMENTA_ACTION_HPP

#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  double momenta_action(eo_ptr<quad_su3> H);
  double momenta_action(quad_su3 *H);
}

#endif
