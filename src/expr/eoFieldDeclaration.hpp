#ifndef _EOFIELDDECLARATION_HPP
#define _EOFIELDDECLARATION_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/old_field.hpp>
#include <base/memory_manager.hpp>
#include <expr/field.hpp>

namespace nissa
{
  /// Structure to hold an even/old field
  template <typename C,
	    typename Fund,
	    FieldLayout FL=defaultFieldLayout,
	    MemoryType MT=defaultMemoryType,
	    bool IsRef=false>
  struct EoField;
}

#endif
