#ifndef _CG_INVERT_STD2EE_M2_HPP
#define _CG_INVERT_STD2EE_M2_HPP

#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void inv_stD2ee_m2_cg(color *sol,color *guess,eo_ptr<quad_su3> conf,double m2,int niter,double residue,color *source);
}

#endif
