#include "new_types/new_types_definitions.hpp"

#ifndef _SOURCE_H
#define _SOURCE_H

namespace nissa
{
  void select_propagator_timeslice(colorspinspin *prop_out,colorspinspin *prop_in,int timeslice);
  void select_propagator_timeslice(su3spinspin *prop_out,su3spinspin *prop_in,int timeslice);
}

#endif
