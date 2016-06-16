#ifndef _CG_INVERT_TMCLOVD_EOPREC_HPP
#define _CG_INVERT_TMCLOVD_EOPREC_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmclovD_cg_eoprec(spincolor *solution_lx,spincolor *guess_Koo,quad_su3 *conf_lx,double kappa,double mass,clover_term_t *Cl_lx,inv_clover_term_t *ext_invCl_lx,int nitermax,double residue,spincolor *source_lx);
}

#endif
