#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "dirac_operator_tmclovQ.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //Apply the Q+ and Q- operator to a spincolor,so that we have Q-^-1 (r==0) and Q+^-1 (r==1) as output
  void reconstruct_tmclov_doublet(LxField<spincolor>& outminus,
				  LxField<spincolor>& outplus,
				  const LxField<quad_su3>& conf,
				  const double& kappa,
				  const LxField<clover_term_t>& Cl,
				  const double& mu,
				  const LxField<spincolor>& in)
  {
    apply_tmclovQ(outminus,conf,kappa,Cl,mu,in);
    
    PAR(0,locVol,
	CAPTURE(mu,
		TO_WRITE(outplus),
		TO_READ(outminus),
		TO_READ(in)),ivol,
      {
	spincolor_copy(outplus[ivol],outminus[ivol]);
	spincolor_summ_the_prod_idouble(outplus[ivol],in[ivol],-2*mu);
      });
  }
}
