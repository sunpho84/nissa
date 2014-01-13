#ifndef GAUGE_SWEEPER_H
#define GAUGE_SWEEPER_H

#include "new_types/new_types_definitions.hpp"
#include "communicate/all_to_all.hpp"

namespace nissa
{
  void init_tlSym_sweeper();
  void init_Wilson_sweeper();
  void init_sweeper(gauge_action_name_t);
}

#endif
