#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "communicate/communicate.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_mix.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "operations/su3_paths/squared_staples.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //Compute the gluonic force for the Wilson plaquette action and summ to the output
  //Passed conf must NOT contain the backfield.
  //Of the result still need to be taken the TA
  THREADABLE_FUNCTION_3ARG(Wilson_force_eo_conf, quad_su3**,F, quad_su3**,eo_conf, double,beta)
  {
    GET_THREAD_ID();
    
    verbosity_lv1_master_printf("Computing Wilson force\n");
    
    double r=beta/3;
    compute_summed_squared_staples_eo_conf(F,eo_conf);
    
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	  for(int mu=0;mu<4;mu++)
	    safe_su3_hermitian_prod_double(F[par][ivol][mu],F[par][ivol][mu],r);
	set_borders_invalid(F[par]);
      }
  }}

  //lx version
  THREADABLE_FUNCTION_3ARG(Wilson_force_lx_conf, quad_su3*,out, quad_su3*,conf, double,beta)
  {
    GET_THREAD_ID();
    
    //compute the squared staples
    compute_summed_squared_staples_lx_conf(out,conf);
    
    //take hermitian*r
    double r=beta/3;
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<4;mu++)
	safe_su3_hermitian_prod_double(out[ivol][mu],out[ivol][mu],r);
  }}
}
