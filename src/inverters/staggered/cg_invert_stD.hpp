#ifndef _CG_INVERT_STD2_M2_H
#define _CG_INVERT_STD2_M2_H

namespace nissa
{
  void inv_stD_cg(color **sol,color *guess,quad_su3 **conf,double m,int niter,double residue,color **source);
  inline void inv_stD_cg(color **sol,quad_su3 **conf,double m,int niter,double residue,color **source)
  {inv_stD_cg(sol,NULL,conf,m,niter,residue,source);}
}

#endif
