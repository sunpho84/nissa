#ifndef _TWISTED_PROPAGATOR_H
#define _TWISTED_PROPAGATOR_H

#include "../../src/new_types/new_types_definitions.h"

void compute_mom_space_twisted_propagator(spinspin *prop,double kappa,double mu,momentum_t theta);
void compute_x_space_twisted_propagator_by_fft(spinspin *prop,double kappa,double mu,momentum_t theta);
void compute_x_space_twisted_propagator_by_inverting(spinspin *prop,double kappa,double mu,momentum_t theta);

#endif
