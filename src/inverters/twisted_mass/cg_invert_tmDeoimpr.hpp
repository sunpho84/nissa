#ifndef _CG_INVERT_TMDEOIMPR_H
#define _CG_INVERT_TMDEOIMPR_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void inv_tmD_cg_eoprec_eos(spincolor *solution_lx,spincolor *guess_Koo,quad_su3 *conf_lx,double kappa,double mu,int nitermax,int rniter,double residue,spincolor *source_lx);
  void inv_tmDkern_eoprec_square_eos(spincolor *sol,spincolor *guess,quad_su3 **conf,double kappa,double mu,int nitermax,int rniter,double residue,spincolor *source);
}

#endif
