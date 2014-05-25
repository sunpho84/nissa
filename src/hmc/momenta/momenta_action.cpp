#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"

namespace nissa
{
  //compute the action of the momenta
  double momenta_action(quad_su3 **H)
  {
    //summ the square of H
    double glb_action_eo[2];
    for(int eo=0;eo<2;eo++)
      double_vector_glb_scalar_prod(&(glb_action_eo[eo]),(double*)(H[eo]),(double*)(H[eo]),4*18*loc_volh);
    
    return (glb_action_eo[EVN]+glb_action_eo[ODD])/2;
  }

  //compute the action of the momenta associated to b
  double B_momenta_action(double *H_B)
  {return (*H_B)*(*H_B)/2;}
}
