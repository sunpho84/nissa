#ifndef _TUPLEAPPLY_HPP
#define _TUPLEAPPLY_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file tupleApply.hpp

#include <tuple>
#include <cstddef>

#include <metaprogramming/inline.hpp>

namespace nissa
{
  namespace impl
  {
    template <typename F,
	      typename T,
	      size_t...Is>
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    decltype(auto) _tupleApply(F&& f,
			       T&& tp,
			       std::index_sequence<Is...>)
    {
      return f(std::get<Is>(std::forward<T>(tp))...);
    }
  }
  
  template <typename F,
	    typename T>
  HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
  decltype(auto) tupleApply(F&& f,
			    T&& tp)
  {
    return impl::_tupleApply(std::forward<F>(f),
			     std::forward<T>(tp),
			     std::make_index_sequence<std::tuple_size_v<std::decay_t<T>>>());
  }
}

#endif
