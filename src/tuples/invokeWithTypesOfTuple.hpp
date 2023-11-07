#ifndef _INVOKEWITHTYPESOFTUPLE_HPP
#define _INVOKEWITHTYPESOFTUPLE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file tuples/invokeWithTypesofTuple.hpp

#include <tuple>

#include <metaprogramming/inline.hpp>

namespace nissa
{
  namespace impl
  {
    /// Call a given function passing explicitly all the types of tuple, and passed args
    ///
    /// Forward declaration
    template <typename TP>
    struct _InvokeWithTypesOfTuple;
    
    /// Call a given function passing explicitly all the types of tuple, and passed args
    ///
    /// Internal implementation
    template <typename...TP>
    struct _InvokeWithTypesOfTuple<std::tuple<TP...>>
    {
      /// Call the callable f explicitly passing the tuple type as template parameters, and the passed arguments
      template <typename F,
		typename...Args>
      INLINE_FUNCTION static constexpr CUDA_HOST_AND_DEVICE
      decltype(auto) exec(const F& f,
			  Args&&...args)
      {
	return f.template operator()<TP...>(std::forward<Args>(args)...);
      }
    };
  }
  
  /// Call a given function passing explicitly all the types of tuple, and passed args
  template <typename TP,
	    typename F,
	    typename...Args>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  decltype(auto) invokeWithTypesOfTuple(const F& f,
					Args&&...args)
  {
    return impl::_InvokeWithTypesOfTuple<TP>::exec(f,std::forward<Args>(args)...);
  }
}

#endif
