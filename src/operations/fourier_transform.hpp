#ifndef _FOURIER_TRANSFORM_HPP
#define _FOURIER_TRANSFORM_HPP

namespace nissa
{
  void pass_spinspin_from_mom_to_x_space(spinspin *out,spinspin *in,double *bc);
  void pass_spinspin_from_x_to_mom_space(spinspin *out,spinspin *in,double *bc);
  void pass_spin1prop_from_mom_to_x_space(spin1prop *out,spin1prop *in,double *bc);
  void pass_spin1prop_from_x_to_mom_space(spin1prop *out,spin1prop *in,double *bc);
  void pass_spin1field_from_mom_to_x_space(spin1field *out,spin1field *in,double *bc,bool bar=false);
  void pass_spin1field_from_x_to_mom_space(spin1field *out,spin1field *in,double *bc,bool bar=false);
  void pass_spin_from_mom_to_x_space(spin *out,spin *in,double *bc,bool bar=false);
  void pass_spin_from_x_to_mom_space(spin *out,spin *in,double *bc,bool bar=false);
}

#endif
