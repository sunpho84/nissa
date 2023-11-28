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
  template <TriviallyCopyable Fund>
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
    
    //ScalarWrapFunctor(const ScalarWrapFunctor& oth) =default;
    
    /// Can run on both GPU and CPU as it is trivially copyable
    static constexpr ExecSpace execSpace=
		execOnCPUAndGPU;
    
    /// Evaluate
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    const Fund& operator()() const
    {
      return val;
    }
  };
  
  /// Scalar quantity
  template <TriviallyCopyable F>
  using Scalar=
    FuncExpr<ScalarWrapFunctor<F>,OfComps<>,F>;
  
  /// Creates a Scalar of type F
  template <TriviallyCopyable F>
  constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  Scalar<F> scalar(const F& f)
  {
    return {ScalarWrapFunctor<F>(f),std::make_tuple()};
  }
}

#endif
