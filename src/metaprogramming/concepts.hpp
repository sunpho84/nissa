#ifndef _CONCEPTS_HPP
#define _CONCEPTS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

/// \file metaprogramming/concepts.hpp

#include <tuple>
#include <type_traits>

#include <metaprogramming/inline.hpp>

namespace nissa
{
  /// Concept to catch trivial and standard
  template <typename T>
  concept TrivialAndStandard=
  std::conjunction_v<std::is_trivial<T>,std::is_standard_layout<T>>;
  
  /// Concept to catch trivially copyable for tuple
  template <typename T>
  INLINE_FUNCTION constexpr bool _isTriviallyCopyableTuple(T)
  {
    return false;
  }
  
  /// Concept to catch trivially copyable for tuple
  template <typename...T>
  INLINE_FUNCTION constexpr bool _isTriviallyCopyableTuple(std::tuple<T...>*)
  {
    return std::conjunction_v<std::is_trivially_copyable<T>...>;
  }
  
  /// Predicate to determine if a type is trivally copyavle as tuple
  template <typename T>
  struct _IsTriviallyCopyableTuple
  {
    static constexpr bool value=_isTriviallyCopyableTuple((T*)nullptr);
  };
  
  /// Concept to catch trivially copyable
  template <typename T>
  concept TriviallyCopyable=
  std::disjunction_v<std::is_trivially_copyable<T>,_IsTriviallyCopyableTuple<std::decay_t<T>>>;
  
  /// Concept to catch castability
  template<typename From,
	   typename To>
  concept HasCastOperatorTo=
  requires(From&& a)
  {
    {a.operator To()};
  };
}

#endif
