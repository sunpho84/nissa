#ifndef _SHIFT_HPP
#define _SHIFT_HPP

#include "../../../../src/new_types/new_types_definitions.hpp"

#include "../types/types.hpp"

void unsafe_shift_spin_up(spin *out,spin *in,momentum_t bc,int mu);
void unsafe_shift_spin_dw(spin *out,spin *in,momentum_t bc,int mu);
void unsafe_shift_spin_updw(spin *out,spin *in,momentum_t bc,int mu,int ud);
void shift_spin_updw(spin *out,spin *in,momentum_t bc,int mu,int ud);
void shift_spin_up(spin *out,spin *in,momentum_t bc,int mu);
void shift_spin_dw(spin *out,spin *in,momentum_t bc,int mu);
void shift_spin1field_up(spin1field *out,spin1field *in,momentum_t bc,int mu);
void shift_spin1field_dw(spin1field *out,spin1field *in,momentum_t bc,int mu);
void shift_spinspin_source_updw(spinspin *out,spinspin *in,momentum_t ext_bc,int mu,int ud);
void shift_spinspin_source_up(spinspin *out,spinspin *in,momentum_t ext_bc,int mu);
void shift_spinspin_source_dw(spinspin *out,spinspin *in,momentum_t ext_bc,int mu);
void shift_spinspin_sink_up(spinspin *out,spinspin *in,momentum_t ext_bc,int mu);
void shift_spinspin_sink_dw(spinspin *out,spinspin *in,momentum_t ext_bc,int mu);
void shift_spinspin_source(spinspin *out,spinspin *in,momentum_t ext_bc,coords r);
void compute_x_space_propagator_to_sink_from_source(spinspin prop,spinspin *q_prop,momentum_t bc,coords sink,coords source);

#endif
