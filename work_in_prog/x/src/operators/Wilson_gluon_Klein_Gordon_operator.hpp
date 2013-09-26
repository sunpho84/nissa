#ifndef _WILSON_GLUON_KLEIN_GORDON_OPERATOR_H
#define _WILSON_GLUON_KLEIN_GORDON_OPERATOR_H

#include "../types/types.hpp"

void apply_Wilson_gluon_mom_Klein_Gordon_operator(spin1field *out,spin1field *in,gluon_info gl);
void apply_Wilson_gluon_x_Klein_Gordon_operator(spin1field *out,spin1field *in,gluon_info gl);

#endif
