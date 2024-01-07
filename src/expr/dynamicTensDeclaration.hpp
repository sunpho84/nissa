#ifndef _DYNAMICTENSORDECLARATION_HPP
#define _DYNAMICTENSORDECLARATION_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/dynamicTensorDeclaration.hpp

#include <base/memoryManager.hpp>

namespace nissa
{
  /// Tensor
  ///
  /// Forward declaration
  template <typename C,
	    typename Fund,
	    MemoryType MT=defaultMemoryType,
	    bool IsRef=false>
  struct DynamicTens;
  
  PROVIDE_FEATURE(DynamicTens);
}

#endif
