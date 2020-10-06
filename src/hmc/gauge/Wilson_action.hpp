#ifndef _WILSON_ACTION_HPP
#define _WILSON_ACTION_HPP

namespace nissa
{
  void Wilson_action(double *action,quad_su3 *lx_conf,double beta);
  void Wilson_action(double *action,quad_su3 **eo_conf,double beta);
}

#endif
