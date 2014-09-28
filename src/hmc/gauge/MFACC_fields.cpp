#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/random.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/momenta/MFACC.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //generate Fourier acceleration-related fields
  THREADABLE_FUNCTION_1ARG(generate_MFACC_fields, su3*,pi)
  {
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      herm_put_to_gauss(pi[ivol],&(loc_rnd_gen[ivol]),1);
    set_borders_invalid(pi);
  }
  THREADABLE_FUNCTION_END
  
  //compute the action for the Fourier acceleration-related fields
  double MFACC_fields_action(su3 **phi)
  {
    //summ the square of pi
    double glb_action_lx[2];
    for(int id=0;id<2;id++)
      double_vector_glb_scalar_prod(&(glb_action_lx[id]),(double*)(phi[id]),(double*)(phi[id]),18*loc_vol);
    
    return (glb_action_lx[0]+glb_action_lx[1])/2;
  }
  
  //Evolve Fourier acceleration related fields
  THREADABLE_FUNCTION_5ARG(evolve_MFACC_fields, su3**,phi, quad_su3*,conf, double,kappa, su3**,pi, double,dt)
  {
    verbosity_lv2_master_printf("Evolving Fourier fields, dt=%lg\n",dt);
    
    //allocate
    su3 *F=nissa_malloc("temp",loc_vol+bord_vol,su3);
    su3 *temp=nissa_malloc("temp",loc_vol+bord_vol,su3);
    
    for(int id=0;id<2;id++)
      {
        //compute
        apply_MFACC(temp,conf,kappa,pi[id]);
        apply_MFACC(F,conf,kappa,temp);
        
        //evolve
        double_vector_summ_double_vector_prod_double((double*)(phi[id]),(double*)(phi[id]),(double*)F,dt,loc_vol*18);
      }
    
    nissa_free(F);
    nissa_free(temp);
  }
  THREADABLE_FUNCTION_END
  
  //Evolve Fourier acceleration related momenta
  THREADABLE_FUNCTION_3ARG(evolve_MFACC_momenta, su3**,pi, su3**,phi, double,dt)
  {
    verbosity_lv2_master_printf("Evolving Fourier momenta, dt=%lg\n",dt);
    
    for(int id=0;id<2;id++)
      double_vector_summ_double_vector_prod_double((double*)(pi[id]),(double*)(pi[id]),(double*)(phi[id]),-dt,loc_vol*18);
  }
  THREADABLE_FUNCTION_END
}
