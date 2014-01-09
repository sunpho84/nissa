#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "dirac_operator_tmQ.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //Apply the Q+ and Q- operator to a spincolor,so that we have Q-^-1 (r==0) and Q+^-1 (r==1) as output
  THREADABLE_FUNCTION_6ARG(reconstruct_tm_doublet, spincolor*,outminus, spincolor*,outplus, quad_su3*,conf, double,kappa, double,mu, spincolor*,in)
  {
    GET_THREAD_ID();
    
    apply_tmQ(outminus,conf,kappa,mu,in);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      unsafe_spincolor_summ_with_ifactor(outplus[ivol],outminus[ivol],in[ivol],-2*mu);
    
    set_borders_invalid(outminus);
    set_borders_invalid(outplus);
  }
  THREADABLE_FUNCTION_END
}
