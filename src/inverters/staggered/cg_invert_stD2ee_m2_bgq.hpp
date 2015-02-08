#ifndef _CG_INVERT_STD2EE_M2_BGQ_H
#define _CG_INVERT_STD2EE_M2_BGQ_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void inv_stD2ee_m2_cg_bgq(bi_color *sol,bi_color *guess,bi_oct_su3 **conf,double m2,int niter,double residue,bi_color *source);
  void inv_stD2ee_m2_bicgstab_bgq(bi_color *sol,bi_color *guess,bi_oct_su3 **conf,double m2,int niter,double residue,bi_color *source);
}

#endif
