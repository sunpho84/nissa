#ifndef _CG_INVERT_TMCLOVD_EOPREC_HPP
#define _CG_INVERT_TMCLOVD_EOPREC_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmclovD_cg_eoprec(spincolor *solution_lx,spincolor *guess_Koo,quad_su3 *conf_lx,double kappa,clover_term_t *Cl_lx,inv_clover_term_t *ext_invCl_lx,double cSW,double mass,int nitermax,double residue,spincolor *source_lx);
}

#endif
