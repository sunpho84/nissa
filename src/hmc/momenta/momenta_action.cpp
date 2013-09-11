#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "base/global_variables.h"
#include "linalgs/linalgs.h"
#include "new_types/new_types_definitions.h"
#include "new_types/su3.h"
#include "routines/mpi_routines.h"

//compute the action of the momenta
double momenta_action(quad_su3 **H)
{
  //summ the square of H
  double glb_action_eo[2];
  for(int eo=0;eo<2;eo++)
    double_vector_glb_scalar_prod(&(glb_action_eo[eo]),(double*)(H[eo]),(double*)(H[eo]),4*18*loc_volh);

  return (glb_action_eo[EVN]+glb_action_eo[ODD])/2;
}
