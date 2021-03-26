#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "dirac_operator_tmclovQ.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //Apply the Q+ and Q- operator to a spincolor,so that we have Q-^-1 (r==0) and Q+^-1 (r==1) as output
  void reconstruct_tmclov_doublet(spincolor* outminus,spincolor* outplus,quad_su3* conf,double kappa,clover_term_t* Cl,double mu,spincolor* in)
  {
    
    apply_tmclovQ(outminus,conf,kappa,Cl,mu,in);
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	spincolor_copy(outplus[ivol],outminus[ivol]);
	spincolor_summ_the_prod_idouble(outplus[ivol],in[ivol],-2*mu);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(outminus);
    set_borders_invalid(outplus);
  }
}
