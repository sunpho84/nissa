#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/bench.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //evolve the momenta with force
  void evolve_lx_momenta_with_force(quad_su3* H,quad_su3* F,double dt)
  {
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int mu=0;mu<NDIM;mu++)
        for(int ic1=0;ic1<NCOL;ic1++)
          for(int ic2=0;ic2<NCOL;ic2++)
            complex_subt_the_prod_idouble(H[ivol.nastyConvert()][mu][ic1][ic2],F[ivol.nastyConvert()][mu][ic1][ic2],dt);
    NISSA_PARALLEL_LOOP_END;
    
    THREAD_BARRIER();
  }
  
  //evolve the configuration with the momenta
  void evolve_lx_conf_with_momenta(quad_su3* lx_conf,quad_su3* H,double dt)
  {
    
    verbosity_lv2_master_printf("Evolving conf with momenta, dt=%lg\n",dt);
    
    START_TIMING(conf_evolve_time,nconf_evolve);
    
    //evolve
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  su3 t1,t2;
	  su3_prod_double(t1,H[ivol.nastyConvert()][mu],dt);
	  safe_hermitian_exact_i_exponentiate(t2,t1);
          
	  safe_su3_prod_su3(lx_conf[ivol.nastyConvert()][mu],t2,lx_conf[ivol.nastyConvert()][mu]);
	}
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(lx_conf);
    
    STOP_TIMING(conf_evolve_time);
  }
}
