#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "dirac_operator_tmclovQ.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //Apply the Q+ and Q- operator to a spincolor,so that we have Q-^-1 (r==0) and Q+^-1 (r==1) as output
  THREADABLE_FUNCTION_7ARG(reconstruct_tmclov_doublet, spincolor*,outminus, spincolor*,outplus, quad_su3*,conf, double,kappa, clover_term_t*,Cl, double,mu, spincolor*,in)
  {
    GET_THREAD_ID();
    
    apply_tmclovQ(outminus,conf,kappa,Cl,mu,in);
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
