#ifndef _FOURIER_TRANSFORM_HPP
#define _FOURIER_TRANSFORM_HPP

namespace nissa
{
  void pass_spinspin_from_mom_to_x_space(spinspin *out,spinspin *in,momentum_t bc);
  void pass_spinspin_from_x_to_mom_space(spinspin *out,spinspin *in,momentum_t bc);
  void pass_spin1prop_from_mom_to_x_space(spin1prop *out,spin1prop *in,momentum_t bc);
  void pass_spin1prop_from_x_to_mom_space(spin1prop *out,spin1prop *in,momentum_t bc);
  void pass_spin1field_from_mom_to_x_space(spin1field *out,spin1field *in,momentum_t bc,bool bar=false);
  void pass_spin1field_from_x_to_mom_space(spin1field *out,spin1field *in,momentum_t bc,bool bar=false);
  void pass_spin_from_mom_to_x_space(spin *out,spin *in,momentum_t bc,bool bar=false);
  void pass_spin_from_x_to_mom_space(spin *out,spin *in,momentum_t bc,bool bar=false);
}

#endif
