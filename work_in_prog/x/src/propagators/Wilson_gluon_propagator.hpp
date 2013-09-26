#ifndef _WILSON_GLUON_PROPAGATOR_H
#define _WILSON_GLUON_PROPAGATOR_H

#include "../types/types.hpp"

void compute_mom_space_Wilson_gluon_propagator(spin1prop *prop,gluon_info gl);
void compute_x_space_Wilson_gluon_propagator_by_fft(spin1prop *prop,gluon_info gl);
void multiply_by_x_space_Wilson_gluon_propagator_by_inv(spin1prop *prop_out,spin1prop *prop_in,gluon_info gl);
void compute_x_space_Wilson_gluon_propagator_by_inv(spin1prop *prop,gluon_info gl);

#endif
