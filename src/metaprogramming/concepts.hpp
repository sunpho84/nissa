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
  std::is_trivial_v<T> and std::is_standard_layout_v<T>;
  
  /// Concept to catch trivially copyable
  template <typename T>
  concept TriviallyCopyable=
  std::is_trivially_copyable_v<T>;
  
}

#endif
