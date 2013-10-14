#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "cg_invert_evn_stD.hpp"

#include "base/global_variables.hpp"
#include "dirac_operators/dirac_operator_stD/dirac_operator_stD.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"

//this is the famous trick to invert the full D matrix using e/o precond: sol[ODD]=1/m*(source[ODD]-Doe*sol[EVN])

namespace nissa
{
  void inv_stD_cg(color **sol,quad_su3 **conf,double m,int niter,double residue,color **source)
  {
    inv_evn_stD_cg(sol[EVN],conf,m,niter,residue,source);
    apply_st2Doe(sol[ODD],conf,sol[EVN]);
    double_vector_linear_comb((double*)(sol[ODD]),(double*)(source[ODD]),1/m,(double*)(sol[ODD]),-0.5/m,loc_volh*6);
  }
}
