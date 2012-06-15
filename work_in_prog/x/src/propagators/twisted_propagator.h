#ifndef _TWISTED_PROPAGATOR_H
#define _TWISTED_PROPAGATOR_H

#include "../../../../src/new_types/new_types_definitions.h"

#include "../types/types.h"

void mom_space_twisted_propagator_of_imom(spin1prop prop,quark_info qu,int imom);
void multiply_from_left_by_mom_space_twisted_propagator(spin *out,spin *in,quark_info qu);
void compute_mom_space_twisted_propagator(spinspin *prop,quark_info qu);
void compute_x_space_twisted_propagator_by_fft(spinspin *prop,quark_info qu);
void multiply_from_left_by_x_space_twisted_propagator_by_fft(spin *prop,spin *ext_source,quark_info qu);
void multiply_from_left_by_x_space_twisted_propagator_by_fft(spinspin *prop,spinspin *ext_source,quark_info qu);
void multiply_from_left_by_x_space_twisted_propagator_by_inv(spin *prop,spin *ext_source,quark_info qu);
void multiply_from_left_by_x_space_twisted_propagator_by_inv(spinspin *prop,spinspin *ext_source,quark_info qu);
void compute_x_space_twisted_propagator_by_inv(spinspin *prop,quark_info qu);

#endif
