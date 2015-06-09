#ifndef _FOURIER_TRANSFORM_HPP
#define _FOURIER_TRANSFORM_HPP

namespace nissa
{
  void pass_spinspin_from_mom_to_x_space_source_or_sink(spinspin *out,spinspin *in,int *dirs,double *bc,int source_or_sink);
  inline void pass_spinspin_from_mom_to_x_space_source_or_sink(spinspin *out,spinspin *in,double *bc,int source_or_sink)
  {int dirs[4]={1,1,1,1};pass_spinspin_from_mom_to_x_space_source_or_sink(out,in,dirs,bc,source_or_sink);}
  void pass_spinspin_from_x_to_mom_space_source_or_sink(spinspin *out,spinspin *in,int *dirs,double *bc,int source_or_sink);
  inline void pass_spinspin_from_x_to_mom_space_source_or_sink(spinspin *out,spinspin *in,double *bc,int source_or_sink)
  {int dirs[4]={1,1,1,1};pass_spinspin_from_x_to_mom_space_source_or_sink(out,in,dirs,bc,source_or_sink);}
  void pass_spin1prop_from_mom_to_x_space(spin1prop *out,spin1prop *in,int *dirs,double *bc);
  inline void pass_spin1prop_from_mom_to_x_space(spin1prop *out,spin1prop *in,double *bc)
  {int dirs[4]={1,1,1,1};pass_spin1prop_from_mom_to_x_space(out,in,dirs,bc);}
  void pass_spin1prop_from_x_to_mom_space(spin1prop *out,spin1prop *in,int *dirs,double *bc);
  inline void pass_spin1prop_from_x_to_mom_space(spin1prop *out,spin1prop *in,double *bc)
  {int dirs[4]={1,1,1,1};pass_spin1prop_from_x_to_mom_space(out,in,dirs,bc);}
  void pass_spin1field_from_mom_to_x_space(spin1field *out,spin1field *in,int *dirs,double *bc,bool bar=false);
  inline void pass_spin1field_from_mom_to_x_space(spin1field *out,spin1field *in,double *bc,bool bar=false)
  {int dirs[4]={1,1,1,1};pass_spin1field_from_mom_to_x_space(out,in,dirs,bc,bar);}
  void pass_spin1field_from_x_to_mom_space(spin1field *out,spin1field *in,int *dirs,double *bc,bool bar=false);
  inline void pass_spin1field_from_x_to_mom_space(spin1field *out,spin1field *in,double *bc,bool bar=false)
  {int dirs[4]={1,1,1,1};pass_spin1field_from_x_to_mom_space(out,in,dirs,bc,bar);}
  void pass_spin_from_mom_to_x_space_source_or_sink(spin *out,spin *in,int *dirs,double *bc,int source_or_sink);
  inline void pass_spin_from_mom_to_x_space_source_or_sink(spin *out,spin *in,double *bc,int source_or_sink)
  {int dirs[4]={1,1,1,1};pass_spin_from_mom_to_x_space_source_or_sink(out,in,dirs,bc,source_or_sink);}
  void pass_spin_from_x_to_mom_space_source_or_sink(spin *out,spin *in,int *dirs,double *bc,int source_or_sink);
  inline void pass_spin_from_x_to_mom_space_source_or_sink(spin *out,spin *in,double *bc,int source_or_sink)
  {int dirs[4]={1,1,1,1};pass_spin_from_x_to_mom_space_source_or_sink(out,in,dirs,bc,source_or_sink);}
}

#endif
