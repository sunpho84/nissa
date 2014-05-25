#ifndef _MOMENTA_ACTION_H
#define _MOMENTA_ACTION_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  double momenta_action(quad_su3 **H);
  double B_momenta_action(double *H_B);
}

#endif
