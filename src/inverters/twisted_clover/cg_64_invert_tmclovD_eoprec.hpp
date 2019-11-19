#ifndef _CG_64_INVERT_TMCLOVD_EOPREC_HPP
#define _CG_64_INVERT_TMCLOVD_EOPREC_HPP

#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmclovDkern_eoprec_square_eos_cg_64(spincolor *sol,spincolor *guess,eo_ptr<quad_su3> eo_conf,double kappa,clover_term_t *Cl_odd,inv_clover_term_t *invCl_evn,double mu,int niter,double residue,spincolor *source);
}

#endif
