#ifndef _FEATURE_HPP
#define _FEATURE_HPP

#include <metaprogramming/inline.hpp>
#include <type_traits>

namespace nissa
{
  /* Usable to recognize a FEAT */
#define PROVIDE_FEATURE(NAME)				\
  struct NAME ## Feat					\
  {							\
  };							\
  							\
  /*! Determine if the type T is deriving from NAME */	\
  template <typename T>					\
  inline constexpr bool is ## NAME =			\
    std::is_base_of_v<NAME ## Feat,			\
		      std::decay_t<T>>;			\
  							\
  template <typename T>					\
  concept DerivedFrom ## NAME = is ## NAME<T>

}

#endif
