#ifndef _CG_INVERT_TMCLOVD_EOPREC_HPP
#define _CG_INVERT_TMCLOVD_EOPREC_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmclovD_cg_eoprec(spincolor *solution_lx,spincolor *guess_Koo,quad_su3 *conf_lx,double kappa,double mu,int nitermax,double residue,spincolor *source_lx);
  void inv_tmclovDkern_eoprec_square_eos(spincolor *sol,spincolor *guess,quad_su3 **conf,double kappa,double mu,int nitermax,double residue,spincolor *source);
}

#endif
