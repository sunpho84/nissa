#ifndef _MOMENTA_GENERATION_HPP
#define _MOMENTA_GENERATION_HPP

#include "geometry/geometry_eo.hpp"
#include "new_types/rat_approx.hpp"

namespace nissa
{
  void generate_hmc_momenta(eo_ptr<quad_su3> H);
  void generate_hmc_momenta(quad_su3 *H);
}

#endif
