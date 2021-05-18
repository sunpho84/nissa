#ifndef _NNARY_EXPR_HPP
#define _NNARY_EXPR_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file nnaryExpr.hpp

#include <tensor/expr.hpp>
#include <tensor/exprRefOrVal.hpp>

namespace nissa
{
  /// Nnary expression needed to capture n arguments
  template <typename T,   // Derived class
	    typename _Es, // Nested expressions, in a tuple as passed to constructor
	    typename TC,
	    typename F>
  struct NnaryExpr;
  
  /// Nnary expression needed to capture N arguments
  template <typename T,     // Derived class
	    typename..._Es, // Nested expressions, as it is passed to constructor
	    typename TC,
	    typename F>
  struct NnaryExpr<T,std::tuple<_Es...>,TC,F> :
    Expr<T,TC,F>
  {
    /// Type of the nested expressions
    using NestedExprs=
      std::tuple<ExprRefOrVal<_Es>...>;
    
    /// First nested expression
    NestedExprs nestedExprs;
    
#define DECLARE_NESTED_EXPR(ATTRIB)	\
    /*! Returns the I-th nested expression */	\
    template <size_t I>				\
    decltype(auto) nestedExpr() ATTRIB	\
    {						\
      return					\
	std::get<I>(nestedExprs);		\
    }
    
    DECLARE_NESTED_EXPR(const);
    
    DECLARE_NESTED_EXPR(/* non const */);
    
#undef DECLARE_NESTED_EXPR
    
    /// Returns the I-th nested expression type
    template <size_t I>
    using NestedExpr=
      std::tuple_element_t<I,NestedExprs>;
    
    /// Constructor for the nested expression
    template <typename...U>
    NnaryExpr(U&&...u) :
      nestedExprs(std::forward<U>(u)...)
    {
    }
    
    /// Import assignemnt operator
    using Expr<T,TC,F>::operator=;
  };
}

#endif
