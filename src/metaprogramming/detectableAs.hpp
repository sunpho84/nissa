#ifndef _DETECTABLEAS_HPP
#define _DETECTABLEAS_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file metaprogramming/detectableAs.hpp
///
/// \brief While waiting for concepts, we define a detectable pattern

#include <type_traits>

namespace nissa
{
#define PROVIDE_DETECTABLE_AS(NAME)			\
  struct DetectableAs ## NAME				\
  {							\
  };							\
							\
  template <typename T>					\
  constexpr bool is ## NAME=				\
    std::is_base_of_v<DetectableAs ## NAME,		\
    std::decay_t<T>>
}

#endif
