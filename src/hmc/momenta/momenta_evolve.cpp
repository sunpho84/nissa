#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "inverters/momenta/cg_invert_MFACC.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //accelerate the momenta
  void accelerate_lx_momenta(quad_su3 *M,quad_su3 *conf,double kappa,int niter,double residue,quad_su3 *H)
  {inv_MFACC_cg(M,NULL,conf,kappa,niter,residue,H);}
  
  //evolve the momenta with force
  THREADABLE_FUNCTION_3ARG(evolve_lx_momenta_with_force, quad_su3*,H, quad_su3*,F, double,dt)
  {
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<4;mu++)
        for(int ic1=0;ic1<3;ic1++)
          for(int ic2=0;ic2<3;ic2++)
            complex_subt_the_prod_idouble(H[ivol][mu][ic1][ic2],F[ivol][mu][ic1][ic2],dt);
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END

  //evolve the configuration with the momenta
  THREADABLE_FUNCTION_3ARG(evolve_lx_conf_with_momenta, quad_su3*,lx_conf, quad_su3*,H, double,dt)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Evolving conf with momenta, dt=%lg\n",dt);
    
    //evolve
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<4;mu++)
	{
	  su3 t1,t2;
	  su3_prod_double(t1,H[ivol][mu],dt);
	  safe_anti_hermitian_exact_i_exponentiate(t2,t1);
          
	  safe_su3_prod_su3(lx_conf[ivol][mu],t2,lx_conf[ivol][mu]);
	}
    
    set_borders_invalid(lx_conf);
  }
  THREADABLE_FUNCTION_END
  
  //accelerate and evolve
  void evolve_lx_conf_with_accelerated_momenta(quad_su3 *lx_conf,quad_su3 *H,double kappa,int niter,double residue,double dt)
  {
    if(fabs(kappa)>1.e-10)
      {
	quad_su3 *M=nissa_malloc("M",loc_vol+bord_vol,quad_su3);
	accelerate_lx_momenta(M,lx_conf,kappa,niter,residue,H);
	evolve_lx_conf_with_momenta(lx_conf,M,dt);
	nissa_free(M);
      }
    else evolve_lx_conf_with_momenta(lx_conf,H,dt);
  }
}
