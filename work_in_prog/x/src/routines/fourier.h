#ifndef _FOURIER_H
#define _FOURIER_H

#include "../../../../src/new_types/new_types_definitions.h"

void pass_spinspin_from_mom_to_x_space(spinspin *out,spinspin *in,momentum_t bc);
void pass_spinspin_from_x_to_mom_space(spinspin *out,spinspin *in,momentum_t bc);
void pass_spin1prop_from_mom_to_x_space(spin1prop *out,spin1prop *in,momentum_t bc);
void pass_spin1prop_from_x_to_mom_space(spin1prop *out,spin1prop *in,momentum_t bc);
void pass_spin1field_from_mom_to_x_space(spin1field *out,spin1field *in,momentum_t bc);
void pass_spin1field_from_x_to_mom_space(spin1field *out,spin1field *in,momentum_t bc);

#endif
