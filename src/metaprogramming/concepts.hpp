#ifndef _CONCEPTS_HPP
#define _CONCEPTS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

/// \file metaprogramming/concepts.hpp

#include <type_traits>

namespace nissa
{
  /// Concept to catch trivial and standard
  template <typename T>
  concept TrivialAndStandard=
  std::conjunction_v<std::is_trivial<T>,std::is_standard_layout<T>>;
  
  /// Concept to catch trivially copyable
  template <typename T>
  concept TriviallyCopyable=
  std::is_trivially_copyable_v<T>;
  
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
