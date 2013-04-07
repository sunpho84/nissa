#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/communicate.h"
#include "../../base/global_variables.h"
#include "../../new_types/new_types_definitions.h"
#include "../../operations/su3_paths/rectangles.h"
#include "../../routines/ios.h"

//compute the tree level Symanzik action
void tree_level_Symanzik_action(double *action,quad_su3 **conf,double beta,int stagphases_present)
{
  verbosity_lv1_master_printf("Computing tree level Symanzik action\n");
  
  //coefficient of rectangles and squares
  double b1=-1.0/12,b0=1-8*b1;
  
  //compute shapes
  complex glb_shapes;
  global_plaquette_and_rectangles_eo_conf(glb_shapes,conf);
  if(stagphases_present) glb_shapes[RE]*=-1; //stag phases add (-1)^area
  
  //compute the total action
  (*action)=(b0*6*glb_vol*(1-glb_shapes[RE])+b1*12*glb_vol*(1-glb_shapes[IM]))*beta;
}
