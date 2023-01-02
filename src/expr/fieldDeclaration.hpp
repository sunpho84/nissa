#ifndef _FIELDDECLARATION_HPP
#define _FIELDDECLARATION_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/fieldDeclaration.hpp

            #include <base/field.hpp>

#include <base/memory_manager.hpp>
#include <metaprogramming/hasMember.hpp>
#include <tuples/tupleHasType.hpp>

namespace nissa
{
  /// Default presence of halo
  inline constexpr HaloEdgesPresence defaultHaloPresence=
    HaloEdgesPresence::WITHOUT_HALO;
  
  PROVIDE_HAS_MEMBER(fieldLayout);
  
  /// Field, forward declaration
  template <typename InnerComps,
	    typename Fund,
	    SitesCoverage LC=SitesCoverage::FULL_SPACE,
	    FieldLayout FL=defaultFieldLayout,
	    MemoryType MT=defaultMemoryType,
	    bool IsRef=false>
  struct Field2;
}

#endif
