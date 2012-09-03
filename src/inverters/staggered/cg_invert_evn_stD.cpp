#include "../../linalgs/linalgs.h"
#include "../../dirac_operators/dirac_operator_stD/dirac_operator_stD.h"
#include "../../new_types/new_types_definitions.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../inverters/staggered/cg_invert_stD2ee_m2.h"

void inv_evn_stD_cg(color *sol,quad_su3 **conf,double m,int niter,int rniter,double residue,color **source)
{
  color *temp=nissa_malloc("temp",loc_vol+bord_vol,color);
  evn_apply_stD(temp,conf,m,source);
  inv_stD2ee_m2_cg(sol,NULL,conf,m,niter,rniter,residue,temp);
  nissa_free(temp);
}
