#ifndef _CG_INVERT_MFACC_HPP
#define _CG_INVERT_MFACC_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_MFACC_cg(su3* sol,su3* guess,quad_su3* conf,double kappa,int niter,double residue,su3* source);
  void inv_MFACC_cg(quad_su3* sol,quad_su3* guess,quad_su3* conf,double kappa,int niter,double residue,quad_su3* source);
}

#endif
