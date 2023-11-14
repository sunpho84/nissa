#ifndef _FIELDDECLARATION_HPP
#define _FIELDDECLARATION_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/fieldDeclaration.hpp

#include <base/memory_manager.hpp>
#include <metaprogramming/hasMember.hpp>
#include <tuples/tupleHasType.hpp>

namespace nissa
{
  /// Has or not the halo and the edges
  enum class HaloEdgesPresence{WITHOUT_HALO,WITH_HALO,WITH_HALO_EDGES};
  using enum HaloEdgesPresence;
  
  /// Default presence of halo
  inline constexpr HaloEdgesPresence defaultHaloPresence=
    HaloEdgesPresence::WITHOUT_HALO;
  
  /// Memory layout
  enum class FieldLayout{CPU,GPU};
  
  /// Coverage of the field
  enum FieldCoverage{EVEN_SITES,ODD_SITES,FULL_SPACE,EVEN_OR_ODD_SITES};
  
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
	    typename Fund,
	    FieldCoverage FC=FieldCoverage::FULL_SPACE,
	    FieldLayout FL=defaultFieldLayout,
	    MemoryType MT=defaultMemoryType,
	    bool IsRef=false>
  struct Field2;
  
  PROVIDE_FEATURE(Field2);
}

#endif
