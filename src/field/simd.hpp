#ifndef _SIMD_HPP
#define _SIMD_HPP

#include <base/memory_manager.hpp>

namespace nissa
{
  /// Size of the simd registers, in bytes
  constexpr Size simdRegSize=16;
  
  /// Number of elements of a SIMD vector
  template <typename T>
  constexpr Size simdVectNel=simdRegSize/sizeof(T);
}

#endif
