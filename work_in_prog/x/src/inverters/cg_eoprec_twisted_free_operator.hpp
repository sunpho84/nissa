#ifndef _cg_eoprec_twisted_mass_free_operator_H
#define _cg_eoprec_twisted_mass_free_operator_H

#include "../../../../src/new_types/new_types_definitions.hpp"

#include "../types/types.hpp"

void inv_tmDkern_eoprec_square_eos(spin *sol,spin *guess,quark_info qu,int nitermax,double residue,spin *source);
void inv_tmD_cg_eoprec_eos(spin *solution_lx,spin *guess_Koo,quark_info qu,int nitermax,double residue,spin *source_lx);

#endif
