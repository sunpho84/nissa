#ifndef _WILSON_GLUON_PROPAGATOR_H
#define _WILSON_GLUON_PROPAGATOR_H

#include "../types/types.h"

void compute_mom_space_Wilson_gluon_propagator(spin1prop *prop,gluon_info gl);
void compute_x_space_Wilson_gluon_propagator_by_fft(spin1prop *prop,gluon_info gl);
void compute_x_space_Wilson_gluon_propagator_by_inv(spin1prop *prop,gluon_info gl);

#endif
