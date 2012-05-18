#ifndef _FOURIER_H
#define _FOURIER_H

void pass_spinspin_from_mom_to_x_space(spinspin *out,spinspin *in,momentum_t bc);
void pass_spinspin_from_x_to_mom_space(spinspin *out,spinspin *in,momentum_t bc);
void pass_spin1prop_from_mom_to_x_space(spinspin *out,spinspin *in,momentum_t bc);
void pass_spin1prop_from_x_to_mom_space(spinspin *out,spinspin *in,momentum_t bc);

#endif
