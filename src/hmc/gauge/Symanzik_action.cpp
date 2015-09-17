#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "communicate/communicate.hpp"
#include "base/global_variables.hpp"
#include "new_types/new_types_definitions.hpp"
#include "operations/su3_paths/rectangles.hpp"
#include "routines/ios.hpp"

#include "Symanzik_action.hpp"

namespace nissa
{
  //compute the Symanzik action
  void Symanzik_action(double *action,double *glb_shapes,double beta,double C1,bool stagphases_present)
  {
    verbosity_lv2_master_printf("Computing Symanzik action\n");
    (*action)=(get_C0(C1,stagphases_present)*6*glb_vol*(1-glb_shapes[RE])+C1*12*glb_vol*(1-glb_shapes[IM]))*beta;
  }
  
  //lx wrapper
  void Symanzik_action(double *action,quad_su3 *conf,double beta,double C1,bool stagphases_present)
  {
    //compute shapes
    complex glb_shapes;
    global_plaquette_and_rectangles_lx_conf(glb_shapes,conf);
    Symanzik_action(action,(double*)glb_shapes,beta,C1,stagphases_present);
  }
  
  //eo wrapper
  void Symanzik_action(double *action,quad_su3 **conf,double beta,double C1,bool stagphases_present)
  {
    //compute shapes
    complex glb_shapes;
    global_plaquette_and_rectangles_eo_conf(glb_shapes,conf);
    Symanzik_action(action,(double*)glb_shapes,beta,C1,stagphases_present);
  }
}
