#ifndef _WILSON_FORCE_H
#define _WILSON_FORCE_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void Wilson_force_eo_conf(quad_su3 **F,quad_su3 **eo_conf,double beta,bool phase_pres);
  void Wilson_force_lx_conf(quad_su3 *F,quad_su3 *lx_conf,double beta,bool phase_pres);
}

#endif
