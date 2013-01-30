#include "../../new_types/new_types_definitions.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"

#include "../../dirac_operators/dirac_operator_tmclovQ/dirac_operator_tmclovQ.h"
#include "cg_invert_tmclovQ2.h"

void inv_tmclovQ_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double mu,int niter,int rniter,double residue,spincolor *source)
{
  inv_tmclovQ2_cg(sol,NULL,conf,kappa,csw,Pmunu,mu,niter,rniter,residue,source);
  spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
  vector_copy(temp,sol);
  apply_tmclovQ(sol,conf,kappa,csw,Pmunu,mu,temp);
  nissa_free(temp);
}
