#ifndef _FOURIER_TRANSFORM_HPP
#define _FOURIER_TRANSFORM_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/field.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
#define DEFINE_PASS(TYPE)						\
  void pass_ ## TYPE ## _from_mom_to_x_space(LxField<TYPE>& out,	\
					     const LxField<TYPE>& in,	\
					     const WhichDirs& dirs,	\
					     const Momentum& bc,	\
					     const int& source_or_sink,	\
					     const bool& include_phases); \
									\
  inline void pass_ ## TYPE ## _from_mom_to_x_space(LxField<TYPE>& out,	\
						    const LxField<TYPE>& in, \
						    const Momentum& bc, \
						    const int& source_or_sink, \
						    const bool& include_phases) \
  {									\
    pass_ ## TYPE ## _from_mom_to_x_space(out,in,allDirs,bc,source_or_sink,include_phases); \
  }									\
  									\
  void pass_ ## TYPE ## _from_x_to_mom_space(LxField<TYPE>& out,	\
					     const LxField<TYPE>& in,	\
					     const WhichDirs& dirs,	\
					     const Momentum& bc,	\
					     const int& source_or_sink,	\
					     const bool& include_phases); \
  									\
  inline void pass_ ## TYPE ## _from_x_to_mom_space(LxField<TYPE>& out,	\
						    const LxField<TYPE>& in, \
						    const Momentum& bc, \
						    const int& source_or_sink, \
						    const bool& include_phases) \
  {									\
    pass_ ## TYPE ## _from_x_to_mom_space(out,in,allDirs,bc,source_or_sink,include_phases);\
  }
  
  DEFINE_PASS(spinspin);
  DEFINE_PASS(spin1prop);
  DEFINE_PASS(spin);
  DEFINE_PASS(spincolor);
  DEFINE_PASS(spin1field);
}

#endif
