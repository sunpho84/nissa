#ifndef _DYNAMICTENSORDECLARATION_HPP
#define _DYNAMICTENSORDECLARATION_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/dynamicTensorDeclaration.hpp

#include <base/memory_manager.hpp>
#include <metaprogramming/detectableAs.hpp>

namespace nissa
{
  PROVIDE_DETECTABLE_AS(DynamicTens);
  
  /// Tensor
  ///
  /// Forward declaration
  template <typename C,
	    typename Fund,
	    MemoryType MT=MemoryType::CPU,
	    bool IsRef=false>
  struct DynamicTens;
}

#endif
