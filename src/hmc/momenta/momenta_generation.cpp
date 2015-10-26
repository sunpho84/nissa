#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/random.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "hmc/gauge/MFACC_fields.hpp"
#include "inverters/momenta/cg_invert_MFACC.hpp"
#include "inverters/momenta/cgm_invert_MFACC.hpp"
#include "new_types/rat_approx.hpp"
#include "new_types/su3.hpp"
#include "operations/remez/remez_algorithm.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //generate momenta using guassian hermitean matrix generator
  THREADABLE_FUNCTION_1ARG(generate_hmc_momenta, quad_su3**,H)
  {
    GET_THREAD_ID();
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	  for(int mu=0;mu<NDIM;mu++)
	    herm_put_to_gauss(H[par][ivol][mu],&(loc_rnd_gen[loclx_of_loceo[par][ivol]]),1);
	set_borders_invalid(H[par]);
      }
  }
  THREADABLE_FUNCTION_END
  
  //generate momenta using guassian hermitean matrix generator
  THREADABLE_FUNCTION_3ARG(generate_hmc_momenta, quad_su3*,H, quad_su3*,conf, double,kappa)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++) herm_put_to_gauss(H[ivol][mu],&(loc_rnd_gen[ivol]),1);
    set_borders_invalid(H);
    
    //better idea: generate mu by mu and invert, put in position at the end
    
    //get the rational approx
    rat_approx_t rat_exp_H;
    generate_approx_of_maxerr(rat_exp_H,1e-5,1.0,1e-9,-1,2);
    master_printf_rat_approx(&rat_exp_H);
    //summ_src_and_all_inv_MFACC_cgm(H,conf,kappa,rat_exp_H,1000000,1e-13,in_H);
    
    
    crash("");
  }
  THREADABLE_FUNCTION_END
  
  //generate momenta needed for Fourier acceleration
  THREADABLE_FUNCTION_4ARG(generate_MFACC_momenta, su3**,pi, quad_su3*,conf, double,kappa, double,residue)
  {
    verbosity_lv1_master_printf("Generating Fourier acceleration momenta\n");
    
    //allocate gaussian field
    su3 *V=nissa_malloc("V",loc_vol+bord_vol,su3);
    
    for(int id=0;id<2;id++)
      {
        //generate gaussianly V and then invert on it
        generate_MFACC_fields(V);
        inv_MFACC_cg(pi[id],NULL,conf,kappa,10000000,residue,V);
      }
    verbosity_lv2_master_printf("\n");
    
    nissa_free(V);
  }
  THREADABLE_FUNCTION_END
}
