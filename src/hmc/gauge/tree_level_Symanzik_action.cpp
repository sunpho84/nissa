#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "communicate/communicate.hpp"
#include "base/global_variables.hpp"
#include "new_types/new_types_definitions.hpp"
#include "operations/su3_paths/rectangles.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //compute the tree level Symanzik action
  void tree_level_Symanzik_action(double *action,double *glb_shapes,double beta,bool stagphases_present)
  {
    //coefficient of rectangles and squares
    double tlSym_b1=-1.0/12,tlSym_b0=1-8*tlSym_b1;
  
    verbosity_lv2_master_printf("Computing tree level Symanzik action\n");
    
    //compute the total action
    if(stagphases_present) glb_shapes[RE]*=-1; //stag phases add (-1)^area
    (*action)=(tlSym_b0*6*glb_vol*(1-glb_shapes[RE])+tlSym_b1*12*glb_vol*(1-glb_shapes[IM]))*beta;
  }

  //lx wrapper
  void tree_level_Symanzik_action(double *action,quad_su3 *conf,double beta,bool stagphases_present)
  {
    //compute shapes
    complex glb_shapes;
    global_plaquette_and_rectangles_lx_conf(glb_shapes,conf);
    tree_level_Symanzik_action(action,(double*)glb_shapes,beta,stagphases_present);
  }
  
  //eo wrapper
  void tree_level_Symanzik_action(double *action,quad_su3 **conf,double beta,bool stagphases_present)
  {
    //compute shapes
    complex glb_shapes;
    global_plaquette_and_rectangles_eo_conf(glb_shapes,conf);
    tree_level_Symanzik_action(action,(double*)glb_shapes,beta,stagphases_present);
  }
}
