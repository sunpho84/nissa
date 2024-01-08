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
  // namespace impl
  // {
  //   /// Call a given function passing explicitly all the types of tuple, and passed args
  //   ///
  //   /// Forward declaration
  //   template <typename TP>
  //   struct _InvokeWithTypesOfTuple;
    
  //   /// Call a given function passing explicitly all the types of tuple, and passed args
  //   ///
  //   /// Internal implementation
  //   template <typename...TP>
  //   struct _InvokeWithTypesOfTuple<std::tuple<TP...>>
  //   {
  //     /// Call the callable f explicitly passing the tuple type as template parameters, and the passed arguments
  //     template <typename F,
  // 		typename...Args>
  //     INLINE_FUNCTION static constexpr HOST_DEVICE_ATTRIB
  //     decltype(auto) exec(const F& f,
  // 			  Args&&...args)
  //     {
  // 	return ;
  //     }
  //   };
  // }
  
  /// Call a given function passing explicitly all the types of tuple, and passed args
  template <typename TP,
	    typename F,
	    typename...Args>
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
  decltype(auto) invokeWithTypesOfTuple(const F& f,
					Args&&...args)
  {
    return
      [&f]<typename...T,typename...A>(std::tuple<T...>*,A&&...a) ->decltype(auto)
      {
	return f.template operator()<T...>(std::forward<A>(a)...);
      }((TP*)nullptr,std::forward<Args>(args)...);
  }
}

#endif
