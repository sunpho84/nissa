#ifndef _MOMENTA_ACTION_HPP
#define _MOMENTA_ACTION_HPP

#include "base/field.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  double momenta_action(const EoField<quad_su3>& H);
  
  double momenta_action(quad_su3 *H);
}

#endif
