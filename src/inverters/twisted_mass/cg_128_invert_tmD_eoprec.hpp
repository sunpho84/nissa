#ifndef _CG_128_INVERT_TMDKERN_EOPREC_SQUARE_EOS_HPP
#define _CG_128_INVERT_TMDKERN_EOPREC_SQUARE_EOS_HPP

#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmDkern_eoprec_square_eos_cg_128(spincolor *sol,spincolor *guess,eo_ptr<quad_su3> conf,double kappa,double mu,int niter,double external_solver_residue,spincolor *external_source);
}

#endif
