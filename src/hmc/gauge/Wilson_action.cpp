#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "communicate/communicate.hpp"
#include "base/global_variables.hpp"
#include "new_types/new_types_definitions.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //compute the tree level Symanzik action
  void Wilson_action(double *action,double plaq,double beta,bool stagphases_present)
  {
    verbosity_lv2_master_printf("Computing Wilson gauge action\n");
    if(stagphases_present) plaq*=-1; //stag phases add (-1)^area
    (*action)=beta*6*(1-plaq)*glb_vol;
  }
  
  //eo wrapper
  void Wilson_action(double *action,quad_su3 **eo_conf,double beta,bool stagphases_present)
  {Wilson_action(action,global_plaquette_eo_conf(eo_conf),beta,stagphases_present);}

  //lx wrapper
  void Wilson_action(double *action,quad_su3 *lx_conf,double beta,bool stagphases_present)
  {Wilson_action(action,global_plaquette_lx_conf(lx_conf),beta,stagphases_present);}
}
