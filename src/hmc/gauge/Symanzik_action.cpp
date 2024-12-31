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
  double Symanzik_action(const double* glb_shapes,
			 const double& beta,
			 const double& C1)
  {
    VERBOSITY_LV2_MASTER_PRINTF("Computing Symanzik action\n");
    
    return (get_C0(C1)*6*glbVol*(1-glb_shapes[RE])+C1*12*glbVol*(1-glb_shapes[IM]))*beta;
  }
  
  //lx wrapper
  double Symanzik_action(const LxField<quad_su3>& conf,
			 const double& beta,
			 const double& C1)
  {
    //compute shapes
    complex glb_shapes;
    global_plaquette_and_rectangles_lx_conf(glb_shapes,conf);
    
    return Symanzik_action((double*)glb_shapes,beta,C1);
  }
  
  //eo wrapper
  double Symanzik_action(const EoField<quad_su3>& eo_conf,
			 const double& beta,
			 const double& C1)
  {
    //compute shapes
    complex glb_shapes;
    global_plaquette_and_rectangles_eo_conf(glb_shapes,eo_conf);

    return Symanzik_action((double*)glb_shapes,beta,C1);
  }
}
