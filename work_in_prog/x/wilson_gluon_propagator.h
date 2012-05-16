#ifndef _WILSON_GLUON_PROPAGATOR_H
#define _WILSON_GLUON_PROPAGATOR_H

void compute_mom_space_wilson_gluon_propagator(spin1prop *prop,double alpha,momentum_t bc);
void compute_x_space_wilson_gluon_propagator_by_fft(spin1prop *prop,double alpha,momentum_t bc);

#endif
