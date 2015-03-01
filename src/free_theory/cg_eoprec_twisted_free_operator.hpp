#ifndef _CG_EOPREC_TWISTED_MASS_FREE_OPERATOR_HPP
#define _CG_EOPREC_TWISTED_MASS_FREE_OPERATOR_HPP

#include "new_types/new_types_definitions.hpp"

#include "free_theory_types.hpp"

namespace nissa
{
  void inv_tmDkern_eoprec_square_eos(spin *sol,spin *guess,quark_info qu,int nitermax,double residue,spin *source);
  void inv_tmD_cg_eoprec_eos(spin *solution_lx,spin *guess_Koo,quark_info qu,int nitermax,double residue,spin *source_lx);
}

#endif
