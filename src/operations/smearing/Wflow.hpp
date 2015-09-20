#ifndef _WFLOW_HPP
#define _WFLOW_HPP

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void adaptative_stout_lx_conf(quad_su3 *conf,double *t,double Tmax,double *ext_dt);
  void Wflow_lx_conf(quad_su3 *conf,double dt);
}

#endif
