#ifndef _CG_INVERT_STD2EE_M2_H
#define _CG_INVERT_STD2EE_M2_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void inv_stD2ee_m2_cg(color *sol,color *guess,quad_su3 **conf,double m2,int niter,int rniter,double residue,color *source);
}

#endif
