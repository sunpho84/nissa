#ifndef _CG_INVERT_EVN_STD_HPP
#define _CG_INVERT_EVN_STD_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_evn_stD_cg(color *sol,color *guess,quad_su3 **conf,double m,int niter,double residue,color **source);
  inline void inv_evn_stD_cg(color *sol,quad_su3 **conf,double m,int niter,double residue,color **source)
  {inv_evn_stD_cg(sol,NULL,conf,m,niter,residue,source);}

}

#endif
