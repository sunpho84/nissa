#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "geometry/geometry_lx.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //compute the Wilson action
  void Wilson_action(double *action,double plaq,double beta)
  {
    verbosity_lv2_master_printf("Computing Wilson gauge action\n");
    (*action)=beta*6*(1-plaq)*glb_vol;
  }
  
  //eo wrapper
  void Wilson_action(double *action,quad_su3 **eo_conf,double beta)
  {Wilson_action(action,global_plaquette_eo_conf(eo_conf),beta);}
  
  //lx wrapper
  void Wilson_action(double *action,quad_su3 *lx_conf,double beta)
  {Wilson_action(action,global_plaquette_lx_conf(lx_conf),beta);}
}
