#ifndef _FOURIER_TRANSFORM_HPP
#define _FOURIER_TRANSFORM_HPP

#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
#define DEFINE_PASS(TYPE)						\
  void pass_ ## TYPE ## _from_mom_to_x_space(TYPE *out,TYPE *in,const Coords<bool>& dirs,const Momentum& bc,int source_or_sink,bool include_phases); \
  inline void pass_ ## TYPE ## _from_mom_to_x_space(TYPE *out,TYPE *in,const Momentum& bc,int source_or_sink,bool include_phases)	\
  {pass_ ## TYPE ## _from_mom_to_x_space(out,in,all_dirs,bc,source_or_sink,include_phases);} \
  void pass_ ## TYPE ## _from_x_to_mom_space(TYPE *out,TYPE *in,const Coords<bool>& dirs,const Momentum& bc,int source_or_sink,bool include_phases); \
  inline void pass_ ## TYPE ## _from_x_to_mom_space(TYPE *out,TYPE *in,const Momentum& bc,int source_or_sink,bool include_phases)	\
  {pass_ ## TYPE ## _from_x_to_mom_space(out,in,all_dirs,bc,source_or_sink,include_phases);}
  
  DEFINE_PASS(spinspin);
  DEFINE_PASS(spin1prop);
  DEFINE_PASS(spin);
  DEFINE_PASS(spincolor);
  DEFINE_PASS(spin1field);
}

#endif
