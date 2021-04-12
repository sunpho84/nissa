#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <base/vectors.hpp>
#include <inverters/staggered/cg_invert_stD2ee_m2_portable.hpp>
#include <geometry/geometry_eo.hpp>
#include <geometry/geometry_mix.hpp>
#include <routines/ios.hpp>

namespace nissa
{
  void inv_stD2ee_m2_cg(color *sol,color *guess,eo_ptr<quad_su3> eo_conf,double m2,int niter,double residue,color *source)
  {
    inv_stD2ee_m2_cg_portable(sol,guess,eo_conf,m2,niter,residue,source);
  }
}
