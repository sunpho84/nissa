#ifndef _CG_64_INVERT_TMDEOIMPR_HPP
#define _CG_64_INVERT_TMDEOIMPR_HPP

namespace nissa
{
  void inv_tmDkern_eoprec_square_eos_cg_64(spincolor *sol,spincolor *guess,quad_su3 **conf,double kappa,double mu,int niter,double residue,spincolor *source);
}

#endif
