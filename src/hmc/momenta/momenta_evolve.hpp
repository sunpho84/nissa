#ifndef _MOMENTA_EVOLVE_HPP
#define _MOMENTA_EVOLVE_HPP

#include "base/field.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void evolve_lx_momenta_with_force(LxField<quad_su3>& H,
				    const LxField<quad_su3>& F,
				    const double& dt);
  
  void evolve_lx_conf_with_momenta(quad_su3 *lx_conf,quad_su3 *H,double dt);
}

#endif
