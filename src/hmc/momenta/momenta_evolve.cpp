#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/bench.hpp"
#include "base/field.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //evolve the momenta with force
    void evolve_lx_momenta_with_force(LxField<quad_su3>& H,
				    const LxField<quad_su3>& F,
				    const double& dt)
  {
    PAR(0,
	locVol,
	CAPTURE(TO_WRITE(H),
		TO_READ(F),
		dt),
	ivol,
	{
	  for(int mu=0;mu<NDIM;mu++)
	    for(int ic1=0;ic1<NCOL;ic1++)
	      for(int ic2=0;ic2<NCOL;ic2++)
		complex_subt_the_prod_idouble(H[ivol][mu][ic1][ic2],
					      F[ivol][mu][ic1][ic2],
					      dt);
	});
  }
  
  //evolve the configuration with the momenta
  void evolve_lx_conf_with_momenta(LxField<quad_su3>& lx_conf,
				   const LxField<quad_su3>& H,
				   const double& dt)
  {
    VERBOSITY_LV2_MASTER_PRINTF("Evolving conf with momenta, dt=%lg\n",dt);
    
    START_TIMING(conf_evolve_time,nconf_evolve);
    
    //evolve
    PAR(0,
	locVol,
	CAPTURE(TO_WRITE(lx_conf),
		TO_READ(H),
		dt),
	ivol,
	{
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      su3 t1,t2;
	      su3_prod_double(t1,H[ivol][mu],dt);
	      safe_hermitian_exact_i_exponentiate(t2,t1);
	      
	      safe_su3_prod_su3(lx_conf[ivol][mu],t2,lx_conf[ivol][mu]);
	    }
	});
    
    STOP_TIMING(conf_evolve_time);
  }
}
