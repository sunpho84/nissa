#ifndef GAUGE_SWEEPER_H
#define GAUGE_SWEEPER_H

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/new_types_definitions.hpp"
#include "communicate/all_to_all.hpp"

namespace nissa
{
#ifdef BGQ
  void compute_tlSym_staples_packed_bgq(su3 staples1,su3 staples2,bi_su3 *links);
#endif
  void init_tlSym_sweeper();
  void init_Wilson_sweeper();
  void init_sweeper(gauge_action_name_t);
}

#endif
