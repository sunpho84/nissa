#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "base/thread_macros.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "gauge/Wilson_force.hpp"
#include "gauge/tree_level_Symanzik_force.hpp"
#include "backfield.hpp"

namespace nissa
{
  //Finish the computation multiplying for the conf and taking TA
  THREADABLE_FUNCTION_2ARG(gluonic_force_finish_computation, quad_su3*,F, quad_su3*,conf)
  {
    GET_THREAD_ID();
    
    //remove the staggered phase from the conf, since they are already implemented in the force
    addrem_stagphases_to_lx_conf(conf);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<4;mu++)
	{
	  su3 temp;
	  unsafe_su3_prod_su3(temp,conf[ivol][mu],F[ivol][mu]);
	  unsafe_su3_traceless_anti_hermitian_part(F[ivol][mu],temp);
	}
    
    //readd
    addrem_stagphases_to_lx_conf(conf);
  }}

  //compute only the gauge part
  THREADABLE_FUNCTION_3ARG(compute_gluonic_force_lx_conf, quad_su3*,F, quad_su3*,conf, theory_pars_t*,physics)
  {
    GET_THREAD_ID();
    
#ifdef BENCH
    if(IS_MASTER_THREAD)
      {
	nglu_comp++;
	glu_comp_time-=take_time();
      }
#endif
    
    switch(physics->gauge_action_name)
      {
      case Wilson_action: Wilson_force_lx_conf(F,conf,physics->beta);break;
      case tlSym_action: tree_level_Symanzik_force_lx_conf(F,conf,physics->beta);break;
      default: crash("Unknown action");
      }
    
    //add the stag phases to the force term, to cancel the one entering the force
    addrem_stagphases_to_lx_conf(F);
    
    gluonic_force_finish_computation(F,conf);
    
#ifdef BENCH
    if(IS_MASTER_THREAD) glu_comp_time+=take_time();
#endif
  }}
}
