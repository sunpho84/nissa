#ifndef _CG_64_INVERT_TMCLOV_DEOPREC_BGQ_HPP
#define _CG_64_INVERT_TMCLOV_DEOPREC_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmclovDkern_eoprec_square_eos_cg_64_bgq(vir_spincolor *sol,vir_spincolor *guess,vir_oct_su3 **vir_eo_conf,double kappa,vir_clover_term_t *Cl_odd,vir_inv_clover_term_t *invCl_evn,double mu,int niter,double residue,vir_spincolor *source);
}

#endif
