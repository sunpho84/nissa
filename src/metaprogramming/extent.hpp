#ifndef _EXTENT_HPP
#define _EXTENT_HPP

#include <cstddef>
#include <type_traits>

namespace nissa
{
  namespace impl
  {
    /// Counts the number of extents
    ///
    /// Internal implementation
    template <typename U>
    static constexpr size_t _nOfExtents()
    {
      if constexpr(std::is_array_v<U>)
	return 1+_nOfExtents<std::remove_extent_t<U>>();
      else return 0;
    };
  }
  
  /// Counts the number of extents
  template <typename T>
  static constexpr size_t nOfExtents=
    impl::_nOfExtents<T>();
}

#endif
