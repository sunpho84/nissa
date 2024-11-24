#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/rat_approx.hpp"
#include "new_types/su3_op.hpp"
#include "operations/remez/remez_algorithm.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //generate momenta using guassian hermitean matrix generator
  void generate_hmc_momenta(EoField<quad_su3>& H)
  {
    for(int par=0;par<2;par++)
      {
	PAR(0,locVolh,
	    CAPTURE(par,
		    TO_WRITE(H)),
	    ieo,
	    {
	      for(int mu=0;mu<NDIM;mu++)
		herm_put_to_gauss(H[par][ieo][mu],&(loc_rnd_gen[loclx_of_loceo[par][ieo]]),1);
	    });
      }
  }
  
  //similar for lx
  void generate_hmc_momenta(quad_su3* H)
  {
    crash("Reimplement");
    
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   for(int mu=0;mu<NDIM;mu++)
    // 	herm_put_to_gauss(H[ivol][mu],&(loc_rnd_gen[ivol]),1);
    // NISSA_PARALLEL_LOOP_END;
    
    // set_borders_invalid(H);
  }
}
