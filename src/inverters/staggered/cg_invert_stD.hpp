#ifndef _CG_INVERT_STD2_M2_H
#define _CG_INVERT_STD2_M2_H

namespace nissa
{
  void inv_stD_cg(color **sol,quad_su3 **conf,double m,int niter,int rniter,double residue,color **source);
}

#endif
