#ifndef _STOCHASTIC_PROPAGATORS_H
#define _STOCHASTIC_PROPAGATORS_H

#include "../../../../src/new_types/new_types_definitions.h"

#include "../types/types.h"

void generate_stochastic_tlSym_gluon_propagator(spin1field *phi,spin1field *eta,gluon_info gl);
void generate_stochastic_source_and_tlSym_gluon_propagator(spin1field *phi,spin1field *eta,gluon_info gl);

#endif
