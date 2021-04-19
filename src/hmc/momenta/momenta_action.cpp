#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"

namespace nissa
{
  //compute the action of the momenta
  double momenta_action(eo_ptr<quad_su3> H)
  {
    //summ the square of H
    double glb_action_eo[2];
    FOR_BOTH_PARITIES(par)
      double_vector_glb_scalar_prod(&(glb_action_eo[par.nastyConvert()]),(double*)(H[par]),(double*)(H[par]),sizeof(quad_su3)/sizeof(double)*locVolh.nastyConvert());
    
    return (glb_action_eo[EVN.nastyConvert()]+glb_action_eo[ODD.nastyConvert()])/2;
  }
  
  //lx version
  double momenta_action(quad_su3 *H)
  {
    //summ the square of H
    double glb_action_lx;
    double_vector_glb_scalar_prod(&(glb_action_lx),(double*)(H),(double*)(H),sizeof(quad_su3)/sizeof(double)*locVol.nastyConvert());
    
    return glb_action_lx/2;
  }
}
