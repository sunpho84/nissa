#ifndef _CG_INVERT_STD2EE_M2_PORTABLE_HPP
#define _CG_INVERT_STD2EE_M2_PORTABLE_HPP

#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void inv_stD2ee_m2_cg_portable(color *sol,color *guess,eo_ptr<quad_su3> conf,double m2,int niter,double residue,color *source);
  void inv_evn_stD_cg_portable(color *sol,eo_ptr<quad_su3> conf,double m,int niter,double residue,color **source);
}

#endif
