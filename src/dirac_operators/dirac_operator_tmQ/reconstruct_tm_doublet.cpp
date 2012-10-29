#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "dirac_operator_tmQ.h"

//Apply the Q+ and Q- operator to a spincolor,so that we have Q-^-1 (r==0) and Q+^-1 (r==1) as output
void reconstruct_tm_doublet(spincolor *outminus,spincolor *outplus,quad_su3 *conf,double kappa,double mu,spincolor *in)
{
  apply_tmQ(outminus,conf,kappa,mu,in);
  nissa_loc_vol_loop(ivol)
    unsafe_spincolor_summ_with_ifactor(outplus[ivol],outminus[ivol],in[ivol],-2*mu);
  
  set_borders_invalid(outminus);
  set_borders_invalid(outplus);
}
