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
      else
	return 0;
    };
  }
  
  /// Counts the number of extents
  template <typename T>
  static constexpr size_t nOfExtents=
    impl::_nOfExtents<T>();
  
  namespace impl
  {
    /// Reconstruct the extent, no extent case
    template <typename T,
	      typename Fund>
    struct _DuplicateExtents
    {
      using type=Fund;
    };
    
    /// Reconstruct a single extent and call recursively
    template <int N,
	      typename T,
	      typename Fund>
    struct _DuplicateExtents<T[N],Fund>
    {
      using type=
	typename _DuplicateExtents<T,Fund>::type[N];
    };
  }
  
  /// Reconstruct all the extents, if T=D[N1][N2]...[NN] replace it with Fund[N1][N2]...[NN]
  template <typename T,
	    typename Fund>
  using DuplicateExtents=
    typename impl::_DuplicateExtents<T,Fund>::type;
  
  using I=int[5][3];
}

#endif
