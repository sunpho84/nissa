#ifndef _SCALAR_HPP
#define _SCALAR_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/scalar.hpp

#include <expr/funcExpr.hpp>

namespace nissa
{
  /// Wraps a scalar quantity
  template <typename Fund>
  struct ScalarWrapFunctor
  {
    /// Stored value
    const Fund val;
    
    /// Default constructor
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    explicit ScalarWrapFunctor(const Fund& val) :
      val(val)
    {
    }
    
    /// Evaluate
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    const Fund& operator()() const
    {
      return val;
    }
  };
  
  /// Scalar quantity
  template <typename F>
  using Scalar=
    FuncExpr<ScalarWrapFunctor<F>,OfComps<>,F>;
  
  /// Creates a Scalar of type F
  template <typename F>
  Scalar<F> scalar(const F& f)
  {
    return ScalarWrapFunctor<F>(f);
  }
}

#endif
