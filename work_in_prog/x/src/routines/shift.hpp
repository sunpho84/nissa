#ifndef _SHIFT_H
#define _SHIFT_H

#include "../../../../src/new_types/new_types_definitions.hpp"

#include "../types/types.hpp"

void unsafe_shppift_spin_up(spin *out,spin *in,momentum_t bc,int mu);
void unsafe_shppift_spin_dw(spin *out,spin *in,momentum_t bc,int mu);
void unsafe_shppift_spin_updw(spin *out,spin *in,momentum_t bc,int mu,int ud);
void shppift_spin_updw(spin *out,spin *in,momentum_t bc,int mu,int ud);
void shppift_spin_up(spin *out,spin *in,momentum_t bc,int mu);
void shppift_spin_dw(spin *out,spin *in,momentum_t bc,int mu);
void shppift_spin1field_up(spin1field *out,spin1field *in,momentum_t bc,int mu);
void shppift_spin1field_dw(spin1field *out,spin1field *in,momentum_t bc,int mu);
void shppift_spinspin_source_updw(spinspin *out,spinspin *in,momentum_t ext_bc,int mu,int ud);
void shppift_spinspin_source_up(spinspin *out,spinspin *in,momentum_t ext_bc,int mu);
void shppift_spinspin_source_dw(spinspin *out,spinspin *in,momentum_t ext_bc,int mu);
void shppift_spinspin_sink_up(spinspin *out,spinspin *in,momentum_t ext_bc,int mu);
void shppift_spinspin_sink_dw(spinspin *out,spinspin *in,momentum_t ext_bc,int mu);
void shppift_spinspin_source(spinspin *out,spinspin *in,momentum_t ext_bc,coords r);
void compute_x_space_propagator_to_sink_from_source(spinspin prop,spinspin *q_prop,momentum_t bc,coords sink,coords source);

#endif
