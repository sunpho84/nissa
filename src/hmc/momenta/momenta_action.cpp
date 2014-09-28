#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "dirac_operators/momenta/MFACC.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

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

  //lx version
  double momenta_action(quad_su3 *H)
  {
    //summ the square of H
    double glb_action_lx;
    double_vector_glb_scalar_prod(&(glb_action_lx),(double*)(H),(double*)(H),4*18*loc_vol);
    
    return glb_action_lx/2;
  }

  //compute the action for the Fourier acceleration-related momenta
  THREADABLE_FUNCTION_4ARG(MFACC_momenta_action, double*,tot_action, su3**,pi, quad_su3*,conf, double,kappa)
  {
    //allocate temporary field where to store output
    su3 *V=nissa_malloc("V",loc_vol,su3);
    
    double glb_action_id[2];
    for(int id=0;id<2;id++)
      {
        //apply the kernel
        apply_MFACC(V,conf,kappa,pi[id]);
        double_vector_glb_scalar_prod(&(glb_action_id[id]),(double*)V,(double*)V,18*loc_vol);
      }
    
    (*tot_action)=(glb_action_id[0]+glb_action_id[1])/2;
    nissa_free(V);
  }
  THREADABLE_FUNCTION_END
  
  //compute the action of the momenta associated to b
  double B_momenta_action(double *H_B)
  {return (*H_B)*(*H_B)/2;}
}
