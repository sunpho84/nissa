#ifndef _CG_INVERT_STD2LEB_EE_M2_PORTABLE_HPP
#define _CG_INVERT_STD2LEB_EE_M2_PORTABLE_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_stD2Leb_ee_m2_cg_portable(color *sol,color *guess,oct_su3 **conf,double m2,int niter,double residue,color *source);
}

#endif
