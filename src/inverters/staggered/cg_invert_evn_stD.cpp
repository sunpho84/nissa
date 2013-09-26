#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/dirac_operator_stD/dirac_operator_stD.hpp"
#include "inverters/staggered/cg_invert_stD2ee_m2.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_evn_stD_cg(color *sol,quad_su3 **conf,double m,int niter,int rniter,double residue,color **source)
  {
    //apply the dagger ...
    color *temp=nissa_malloc("temp",loc_volh+bord_volh,color);
    evn_apply_stD_dag(temp,conf,m,source);
    
    //and invert the DD^+
    inv_stD2ee_m2_cg(sol,NULL,conf,m*m,niter,rniter,residue,temp);
    
    nissa_free(temp);
  }
}
