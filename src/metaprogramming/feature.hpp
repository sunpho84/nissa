#ifndef _FEATURE_HPP
#define _FEATURE_HPP

#include <metaprogramming/inline.hpp>
#include <type_traits>

namespace nissa
{
  /* Usable to recognize a FEAT */
#define PROVIDE_FEATURE(NAME)			\
  template <typename F>				\
  struct NAME ## Feat				\
  {						\
    /* Cast to derived type */			\
    F* operator->()				\
    {						\
      return (F*)this;				\
    }						\
						\
    /* Const cast to derived type*/		\
    const F* operator->() const			\
    {						\
      return (const F*)this;			\
    }						\
						\
    /* Cast to derived type*/			\
    F& operator*()				\
    {						\
      return *(F*)this;				\
    }						\
						\
    /* Const cast to derived type*/		\
    const F& operator*() const			\
    {						\
      return *(const F*)this;			\
    }						\
  };						\
  						\
  template <typename T>				\
  inline constexpr bool is ## NAME =		\
    std::is_base_of_v<NAME ## Feat<T>,T>
}

#endif
