#ifndef _FOREACHINDEXOFTUPLE_HPP
#define _FOREACHINDEXOFTUPLE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file tuples/forEachIndexOfTuple.hpp

#include <tuple>
#include <cstddef>

#include <metaprogramming/inline.hpp>

namespace nissa
{
  
  /// Call a given function passing explicitly an inex in turn
  template <typename TP,
	    typename F,
	    typename...Args>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  void forEachIndexOfTuple(const F& f)
  {
    [&f]<size_t...I>(std::index_sequence<I...>*)
      {
	(f.template operator()<I>(),...);
      }((std::make_index_sequence<std::tuple_size_v<TP>>*)nullptr);
  }
}

#endif
