#ifndef _DERIVATIVES_H
#define _DERIVATIVES_H

#include "../../../../src/new_types/new_types_definitions.hpp"

#include "../types/types.hpp"

void spin1field_bw_derivative(spin1field *der,spin1field *in,momentum_t bc,int mu);
void spin1field_fw_derivative(spin1field *der,spin1field *in,momentum_t bc,int mu);

#endif
