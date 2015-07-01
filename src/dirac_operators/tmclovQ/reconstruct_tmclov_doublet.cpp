#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "dirac_operator_tmclovQ.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //Apply the Q+ and Q- operator to a spincolor,so that we have Q-^-1 (r==0) and Q+^-1 (r==1) as output
  THREADABLE_FUNCTION_8ARG(reconstruct_tmclov_doublet, spincolor*,outminus, spincolor*,outplus, quad_su3*,conf, double,kappa, double,cSW, as2t_su3*,Pmunu, double,mu, spincolor*,in)
  {
    GET_THREAD_ID();
    
    apply_tmclovQ(outminus,conf,kappa,cSW,Pmunu,mu,in);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	spincolor_copy(outplus[ivol],outminus[ivol]);
	spincolor_summ_the_prod_idouble(outplus[ivol],in[ivol],-2*mu);
      }
    
    set_borders_invalid(outminus);
    set_borders_invalid(outplus);
  }
  THREADABLE_FUNCTION_END
}
