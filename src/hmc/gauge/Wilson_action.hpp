#ifndef _WILSON_ACTION_HPP
#define _WILSON_ACTION_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "operations/su3_paths/plaquette.hpp"
#include "geometry/geometry_eo.hpp"

namespace nissa
{
  //compute the Wilson action
  inline void Wilson_action(double *action,double plaq,double beta)
  {
    verbosity_lv2_master_printf("Computing Wilson gauge action\n");
    (*action)=beta*6*(1-plaq)*glbVol;
  }
  
  //lx wrapper
  inline void Wilson_action(double *action,quad_su3 *lx_conf,double beta)
  {
    Wilson_action(action,global_plaquette_lx_conf(lx_conf),beta);
  }
  
  template <typename C>
  void Wilson_action_eo_conf(double *action,
			     const C& eo_conf,
			     const double& beta)
  {
    Wilson_action(action,global_plaquette_eo_conf(eo_conf),beta);
  }
}

#endif
