#ifndef _MOMENTA_GENERATION_HPP
#define _MOMENTA_GENERATION_HPP

#include "base/field.hpp"

namespace nissa
{
  void generate_hmc_momenta(EoField<quad_su3>& H);
  void generate_hmc_momenta(quad_su3 *H);
}

#endif
