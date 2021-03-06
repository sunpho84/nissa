#ifndef _CG_128_INVERT_TMCLOVDKERN_EOPREC_SQUARE_EOS_HPP
#define _CG_128_INVERT_TMCLOVDKERN_EOPREC_SQUARE_EOS_HPP

#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmclovDkern_eoprec_square_eos_cg_128(spincolor *sol,spincolor *guess,eo_ptr<quad_su3> conf,double kappa,double cSW,clover_term_t *Cl_odd,inv_clover_term_t *invCl_evn,double mu,int niter,double external_solver_residue,spincolor *external_source);
}

#endif
