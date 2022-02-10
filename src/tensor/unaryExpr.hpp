#ifndef _UNARYEXPR_HPP
#define _UNARYEXPR_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file unaryExpr.hpp

#include <tensor/expr.hpp>
#include <tensor/exprRefOrVal.hpp>

namespace nissa
{
  /// Unary expression needed to capture 1 argument
  template <typename T,  // Derived class
	    typename _E, // Nested expression, as it is passed to constructor
	    typename TC,
	    typename F>
  struct UnaryExpr :
    Expr<T,TC,F>
  {
    /// Type of the nested expression
    using NestedExpr=
      ExprRefOrVal<_E>;
    
    /// Nested expression
    NestedExpr nestedExpr;
    
    /// Nested fundamental type
    using NestedFund=
      typename std::remove_reference_t<_E>::Fund;
    
    /// Nested flags
    static constexpr ExprFlags NestedFlags=
      std::remove_reference_t<_E>::Flags;
    
    /// Dynamic sizes
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) getDynamicSizes() const
    {
      return nestedExpr.getDynamicSizes();
    }
    
    /// Constructor for the nested expression
    template <typename...U>
    UnaryExpr(U&&...u) :
      nestedExpr(std::forward<U>(u)...)
    {
    }
    
    /// Import assignemnt operator
    using Expr<T,TC,F>::operator=;
  };
}

#endif
