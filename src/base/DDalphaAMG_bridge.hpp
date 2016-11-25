#ifndef _DDALPHAAMG_HPP
#define _DDALPHAAMG_HPP

#include "new_types/su3.hpp"

namespace DD
{
  void init_DDalphaAMG();
  void import_gauge_conf(nissa::quad_su3 *conf);
}

#endif
