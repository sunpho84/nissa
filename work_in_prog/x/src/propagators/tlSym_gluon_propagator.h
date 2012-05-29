#ifndef _TLSYM_GLUON_PROPAGATOR_H
#define _TLSYM_GLUON_PROPAGATOR_H

#include "../types/types.h"

void mom_space_tlSym_gluon_propagator_of_imom(spin1prop prop,gluon_info gl,int imom);
void compute_mom_space_tlSym_gluon_propagator(spin1prop *prop,gluon_info gl);
void multiply_mom_space_tlSym_gluon_propagator(spin1field *out,spin1field *in,gluon_info gl);
void multiply_x_space_tlSym_gluon_propagator_by_fft(spin1prop *out,spin1prop *in,gluon_info gl);
void compute_x_space_tlSym_gluon_propagator_by_fft(spin1prop *prop,gluon_info gl);

#endif
