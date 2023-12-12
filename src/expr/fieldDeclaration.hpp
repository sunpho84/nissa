#ifndef _FIELDDECLARATION_HPP
#define _FIELDDECLARATION_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/fieldDeclaration.hpp

#include <base/memoryType.hpp>
#include <expr/comp.hpp>
#include <expr/comps.hpp>
#include <metaprogramming/hasMember.hpp>
#include <tuples/tupleHasType.hpp>

namespace nissa
{
  /// Has or not the halo
  enum class HaloPresence{WITHOUT_HALO,WITH_HALO};
  using enum HaloPresence;
  
  /// Default presence of halo
  inline constexpr HaloPresence defaultHaloPresence=
    HaloPresence::WITHOUT_HALO;
  
  /// Memory layout
  enum class FieldLayout{CPU,GPU};
  
  /// Predefinite memory layout
  constexpr FieldLayout defaultFieldLayout=
	      FieldLayout::
#ifdef USE_CUDA
	      GPU
#else
	      CPU
#endif
	      ;
  
  /////////////////////////////////////////////////////////////////
  PROVIDE_HAS_MEMBER(fieldLayout);
  
  /// Field, forward declaration
  template <typename InnerComps,
	    typename Fund=double,
	    FieldLayout FL=defaultFieldLayout,
	    MemoryType MT=defaultMemoryType,
	    bool IsRef=false>
  struct Field;
  
  PROVIDE_FEATURE(Field);
  
  /////////////////////////////////////////////////////////////////
  
  DECLARE_PARALLELIZABLE_COMP(LocLxSite,int64_t,locLxSite);
  // DECLARE_PARALLELIZABLE_COMP(LocEoSite,int64_t,locEoSite);
  // DECLARE_PARALLELIZABLE_COMP(LocEvnSite,int64_t,locEvnSite);
  // DECLARE_PARALLELIZABLE_COMP(LocOddSite,int64_t,locOddSite);
  
  DECLARE_DYNAMIC_COMP(LocCoord);
  DECLARE_DYNAMIC_COMP(GlbCoord);
  DECLARE_PARALLELIZABLE_COMP(GlbLxSite,int64_t,glbLxSite);
}

#endif
