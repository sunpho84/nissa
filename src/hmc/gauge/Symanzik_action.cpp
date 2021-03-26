#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"
#include "operations/su3_paths/rectangles.hpp"
#include "routines/ios.hpp"

#include "Symanzik_action.hpp"

namespace nissa
{
  //compute the Symanzik action
  void Symanzik_action(double *action,double *glb_shapes,double beta,double C1)
  {
    verbosity_lv2_master_printf("Computing Symanzik action\n");
    (*action)=(get_C0(C1)*6*glbVol*(1-glb_shapes[RE])+C1*12*glbVol*(1-glb_shapes[IM]))*beta;
  }
  
  //lx wrapper
  void Symanzik_action(double *action,quad_su3 *conf,double beta,double C1)
  {
    //compute shapes
    complex glb_shapes;
    global_plaquette_and_rectangles_lx_conf(glb_shapes,conf);
    Symanzik_action(action,(double*)glb_shapes,beta,C1);
  }
  
  //eo wrapper
  void Symanzik_action(double *action,eo_ptr<quad_su3> conf,double beta,double C1)
  {
    //compute shapes
    complex glb_shapes;
    global_plaquette_and_rectangles_eo_conf(glb_shapes,conf);
    Symanzik_action(action,(double*)glb_shapes,beta,C1);
  }
}
