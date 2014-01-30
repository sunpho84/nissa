#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "operations/su3_paths/topological_charge.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif
#include "routines/ios.hpp"

#include "gluonic_force.hpp"

namespace nissa
{
  //compute the topological charge
  THREADABLE_FUNCTION_3ARG(compute_topological_force_lx_conf, quad_su3*,F, quad_su3*,conf, topotential_pars_t*,pars)
  {
    GET_THREAD_ID();
    
    verbosity_lv1_master_printf("Computing topological force\n");
    
    //compute the staples
    topological_staples(F,conf);
    
    //normalize
    double norm=pars->theta/(M_PI*M_PI*128);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<4;mu++)
	safe_su3_hermitian_prod_double(F[ivol][mu],F[ivol][mu],norm);
    set_borders_invalid(F);
  
    //add the stag phases to the force term, to cancel the one entering the force
    addrem_stagphases_to_lx_conf(F);
    
    //take TA
    gluonic_force_finish_computation(F,conf);
  }
  THREADABLE_FUNCTION_END
}
