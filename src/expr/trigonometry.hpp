#ifndef _TRIGONOMETRY_HPP
#define _TRIGONOMETRY_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/trigonometry.hpp0

#include <expr/cWiseCombine.hpp>

namespace nissa
{
  /// Providing overload to std library routines, as they are
  /// obfuscated by the expression overload
#define OVERLOAD_FUN_TYPE(NAME,TYPE,SUFF)		\
  template <typename F>					\
  double NAME(const F& f)				\
    requires(isSafeNumericConversion<F,TYPE>)		\
  {							\
    return NAME ## SUFF((TYPE)f);			\
  }
  
#define PROVIDE_UNARY_FUNCTION(FUN)			\
							\
  OVERLOAD_FUN_TYPE(FUN,double,h);			\
  							\
  OVERLOAD_FUN_TYPE(FUN,float,f);			\
							\
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
  
  PROVIDE_UNARY_FUNCTION(sin);
  PROVIDE_UNARY_FUNCTION(cos);
  PROVIDE_UNARY_FUNCTION(tan);
  PROVIDE_UNARY_FUNCTION(asin);
  PROVIDE_UNARY_FUNCTION(acos);
  PROVIDE_UNARY_FUNCTION(atan);
  
#undef PROVIDE_UNARY_FUNCTION
  
#undef OVERLOAD_FUN_TYPE
}

#endif
