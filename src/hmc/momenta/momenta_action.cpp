#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/field.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  /// Compute the action of the momenta
  double momenta_action(const EoField<quad_su3>& H)
  {
    //summ the square of H
    double action=0;
    for(int eo=0;eo<2;eo++)
      action+=H[eo].norm2();
    
    return action/2;
  }
  
  //lx version
  double momenta_action(quad_su3 *H)
  {
    CRASH("Reimplement");
    // //summ the square of H
    // double glb_action_lx;
    // double_vector_glb_scalar_prod(&(glb_action_lx),(double*)(H),(double*)(H),sizeof(quad_su3)/sizeof(double)*locVol);
    
    // return glb_action_lx/2;
    return 0;
  }
}
