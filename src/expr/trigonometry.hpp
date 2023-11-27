#ifndef _TRIGONOMETRY_HPP
#define _TRIGONOMETRY_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/trigonometry.hpp0

#include <expr/cWiseCombine.hpp>

namespace nissa
{
#define PROVIDE_UNARY_FUNCTION(FUN)			\
  /*! Calculate the unary expression FUN */		\
  struct FUN ## Functor					\
  {							\
    template <typename E>				\
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE	\
    static auto compute(E&& e)				\
    {							\
      return FUN(std::forward<E>(e));			\
    }							\
  };							\
  							\
  /*! Catch FUN(node) */				\
  template <DerivedFromNode E>				\
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE	\
  auto FUN(E&& e)					\
  {							\
    return						\
      cWiseCombine<FUN ## Functor>(std::forward<E>(e));	\
  }
  
  // PROVIDE_UNARY_FUNCTION(sin);
  // PROVIDE_UNARY_FUNCTION(tan);
  // PROVIDE_UNARY_FUNCTION(asin);
  // PROVIDE_UNARY_FUNCTION(acos);
  // PROVIDE_UNARY_FUNCTION(atan);
  
#undef PROVIDE_UNARY_FUNCTION
}

#endif
