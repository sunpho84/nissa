#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "random/randomGenerate.hpp"
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
  void generate_hmc_momenta(eo_ptr<quad_su3> H)
  {
    FOR_BOTH_PARITIES(par)
      {
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  FOR_ALL_DIRS(mu)
	  herm_put_to_gauss(H[par][ieo.nastyConvert()][mu.nastyConvert()],&(loc_rnd_gen[loclx_of_loceo(par,ieo).nastyConvert()]),1);
	NISSA_PARALLEL_LOOP_END;
	
	set_borders_invalid(H[par]);
      }
  }
  //similar for lx
  void generate_hmc_momenta(quad_su3* H)
  {
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      FOR_ALL_DIRS(mu)
      herm_put_to_gauss(H[ivol.nastyConvert()][mu.nastyConvert()],&(loc_rnd_gen[ivol.nastyConvert()]),1);
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(H);
  }
}
