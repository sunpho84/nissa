#ifndef _ARR_HPP
#define _ARR_HPP

#include <metaprogramming/inline.hpp>

namespace nissa
{
  /// Fixed-size aggregate array with host and device element access
  template <typename T,
            int N>
  struct Arr
  {
    /// Contiguous element storage
    T data[N];
    
    /// Access a constant element
    CUDA_HOST_AND_DEVICE
    const T& operator[](const int i) const
    {
      return data[i];
    }
    
    /// Access a mutable element
    CUDA_HOST_AND_DEVICE
    T& operator[](const int i)
    {
      return data[i];
    }
  };
}

#endif
