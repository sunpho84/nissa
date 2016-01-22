#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "geometry/geometry_lx.hpp"
#include "inverters/staggered/cg_invert_stD2ee_m2.hpp"
#include "linalgs/linalgs.hpp"

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_evn_stD_cg(color *sol,color *guess,quad_su3 **conf,double m,int niter,double residue,color **source)
  {
    //apply the dagger ...
    color *temp=nissa_malloc("temp",loc_volh+bord_volh,color);
    evn_apply_stD_dag(temp,conf,m,source);
    
    //and invert the DD^+
    inv_stD2ee_m2_cg(sol,guess,conf,m*m,niter,residue,temp);
    
    nissa_free(temp);
  }
}
