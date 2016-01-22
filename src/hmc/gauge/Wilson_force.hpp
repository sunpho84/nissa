#ifndef _WILSON_FORCE_HPP
#define _WILSON_FORCE_HPP

namespace nissa
{
  void Wilson_force_eo_conf(quad_su3 **F,quad_su3 **eo_conf,double beta);
  void Wilson_force_lx_conf(quad_su3 *F,quad_su3 *lx_conf,double beta);
}

#endif
