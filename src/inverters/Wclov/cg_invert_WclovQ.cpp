#include "../../new_types/new_types_definitions.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"

#include "../../dirac_operators/dirac_operator_WclovQ/dirac_operator_WclovQ.h"
#include "cg_invert_WclovQ2.h"

void inv_WclovQ_cg(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,int niter,int rniter,double residue,spincolor *source)
{
  inv_WclovQ2_cg(sol,NULL,conf,kappa,csw,Pmunu,niter,rniter,residue,source);
  spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
  vector_copy(temp,sol);
  apply_WclovQ(sol,conf,kappa,csw,Pmunu,temp);
}
