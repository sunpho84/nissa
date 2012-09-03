#include "cg_invert_stD2ee_m2.h"
#include "../../dirac_operators/dirac_operator_stD/dirac_operator_stD.h"
#include "../../new_types/new_types_definitions.h"
#include "../../base/global_variables.h"
#include "../../linalgs/linalgs.h"

//this is the famous trick to invert the full D matrix using e/o precond: sol[ODD]=1/m*(source[ODD]-0.5*Doe*sol[EVN])

void inv_stD_cg(color **sol,quad_su3 **conf,double m,int niter,int rniter,double residue,color **source)
{
  inv_evn_stD_cg(sol[EVN],conf,m,niter,rniter,residue,source);
  apply_st2Doe(sol[ODD],conf,sol[EVN]);
  double_vector_linear_comb((double*)(sol[ODD]),(double*)(source[ODD]),1/m,(double*)(sol[ODD]),-0.5/m,loc_vol*6);
}
