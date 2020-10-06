#ifndef _CG_INVERT_STD2EE_M2_BGQ_HPP
#define _CG_INVERT_STD2EE_M2_BGQ_HPP

#include "new_types/rat_approx.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void inv_stD2ee_m2_cg_bgq(vir_color *sol,vir_color *guess,vir_oct_su3 **conf,double m2,int niter,double residue,vir_color *source);
  void inv_stD2ee_m2_bicgstab_bgq(vir_color *sol,vir_color *guess,vir_oct_su3 **conf,double m2,int niter,double residue,vir_color *source);
}

#endif
