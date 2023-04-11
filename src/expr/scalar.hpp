#ifndef _SCALAR_HPP
#define _SCALAR_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/nodes/scalar.hpp

#include <expr/comps.hpp>
#include <expr/dynamicCompsProvider.hpp>
#include <expr/nodeDeclaration.hpp>

namespace nissa
{
  PROVIDE_DETECTABLE_AS(Scalar);
  
  /// Scalar
  ///
  /// Forward declaration
  template <typename F>
  struct Scalar:
    DynamicCompsProvider<CompsList<>>,
    DetectableAsScalar,
    Node<Scalar<F>>
  {
    using Base=Node<Scalar<F>>;
    
    /// Import base class assigners
    using Base::operator=;
    
    /// Components
    using Comps=CompsList<>;
    
    /// Fundamental type
    using Fund=F;
    
    /// Executes where allocated
    // static constexpr auto execSpace=
    //   ExecSpace::HOST_DEVICE;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    CompsList<> getDynamicSizes() const
    {
      return {};
    }
    
    /// No dynamic component
    using DynamicComps=CompsList<>;
    
    /// Value
    const Fund value;
    
    /// We can easily copy
    static constexpr bool storeByRef=false;
    
    /// Evaluate returning always the value
    template <typename...C>
    constexpr INLINE_FUNCTION
    Fund eval(const C&...) const
    {
      return value;
    }
    
    /// States whether the scalar can be simdified
    static constexpr bool canSimdify=
      false;
    
    /// Components on which simdifying
    using SimdifyingComp=
      void;
    
    /// Gets a copy
    Scalar<F> getRef() const
    {
      return value;
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Returns from another F
    INLINE_FUNCTION
    Scalar<F> recreateFromExprs(F& val) const
    {
      return scalar(val);
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Constructor
    constexpr Scalar(const F& value) :
      value(value)
    {
    }
    
    /// Cast to fundamental
    operator const Fund&() const
    {
      return value;
    }
  };
  
  /// Creates a scalar
  template <typename F>
  constexpr INLINE_FUNCTION
  Scalar<F> scalar(const F& val)
  {
    return Scalar(val);
  }
}

#endif
