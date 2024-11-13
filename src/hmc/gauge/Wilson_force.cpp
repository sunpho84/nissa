#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/field.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "operations/su3_paths/squared_staples.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //Compute the gluonic force for the Wilson plaquette action and summ to the output
  //Passed conf must NOT contain the backfield.
  //Of the result still need to be taken the TA
  void Wilson_force_eo_conf(eo_ptr<quad_su3> F,eo_ptr<quad_su3> eo_conf,double beta)
  {
    crash("reimplement");
    
    // verbosity_lv1_master_printf("Computing Wilson force (eo)\n");
    
    // double r=-beta/NCOL;
    // compute_summed_squared_staples_eo_conf(F,eo_conf);
    
    // for(int par=0;par<2;par++)
    //   {
    // 	PAR(0,locVolh,
    // 	    CAPTURE(TO_READ(F)),
    // 	    ivol,
    // 	  {
    // 	    for(int mu=0;mu<NDIM;mu++)
    // 	      safe_su3_hermitian_prod_double(F[par][ivol][mu],F[par][ivol][mu],r);
    // 	  });
    //   }
  }
  
  //lx version
  void Wilson_force_lx_conf(LxField<quad_su3>& F,
			    const LxField<quad_su3>& conf,
			    const double beta)
  {
    crash("reimplement");
    // verbosity_lv1_master_printf("Computing Wilson force (lx)\n");
    
    // const double r=-beta/NCOL;
    // compute_summed_squared_staples_lx_conf(F,conf);
    
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   {
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  safe_su3_hermitian_prod_double(F[ivol][mu],F[ivol][mu],r);
    //   }
    // NISSA_PARALLEL_LOOP_END;
    
    // set_borders_invalid(F);
  }
}
