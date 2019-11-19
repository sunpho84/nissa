#ifndef _CG_INVERT_STD2_M2_H
#define _CG_INVERT_STD2_M2_H

#include "geometry/geometry_eo.hpp"

namespace nissa
{
  void inv_stD_cg(eo_ptr<color> sol,color *guess,eo_ptr<quad_su3> conf,double m,int niter,double residue,eo_ptr<color> source);
  inline void inv_stD_cg(eo_ptr<color> sol,eo_ptr<quad_su3> conf,double m,int niter,double residue,eo_ptr<color> source)
  {inv_stD_cg(sol,NULL,conf,m,niter,residue,source);}
}

#endif
