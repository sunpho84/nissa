#ifndef _CG_INVERT_EVN_STD_H
#define _CG_INVERT_EVN_STD_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void inv_evn_stD_cg(color *sol,color *guess,quad_su3 **conf,double m,int niter,double residue,color **source);
  inline void inv_evn_stD_cg(color *sol,quad_su3 **conf,double m,int niter,double residue,color **source)
  {inv_evn_stD_cg(sol,NULL,conf,m,niter,residue,source);}

}

#endif
