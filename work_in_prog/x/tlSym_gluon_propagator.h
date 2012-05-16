#ifndef _TLSYM_GLUON_PROPAGATOR_H
#define _TLSYM_GLUON_PROPAGATOR_H

void compute_mom_space_tlSym_gluon_propagator(spin1prop *prop,double alpha,double,momentum_t bc);
void compute_x_space_tlSym_gluon_propagator_by_fft(spin1prop *prop,double alpha,double c1,momentum_t bc);

#endif
