#ifndef _WILSON_ACTION_HPP
#define _WILSON_ACTION_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "operations/su3_paths/plaquette.hpp"
#include "geometry/geometry_eo.hpp"

namespace nissa
{
  /// Compute the Wilson action
  inline double Wilson_action(const double& plaq,
			      const double& beta)
  {
    verbosity_lv2_master_printf("Computing Wilson gauge action\n");
    
    return beta*6*(1-plaq)*glbVol;
  }
  
  //lx wrapper
  inline double Wilson_action(const LxField<quad_su3>&lx_conf,
			      const double& beta)
  {
    return Wilson_action(global_plaquette_lx_conf(lx_conf),beta);
  }
  
  template <typename C>
  double Wilson_action(const C& eo_conf,
		       const double& beta)
  {
    return Wilson_action(global_plaquette_eo_conf(eo_conf),beta);
  }
}

#endif
