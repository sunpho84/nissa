#ifndef _NNARY_EXPR_HPP
#define _NNARY_EXPR_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file nnaryExpr.hpp

#include <base/debug.hpp>
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
    
#define DECLARE_NESTED_EXPR(ATTRIB)		\
    /*! Returns the I-th nested expression */	\
    template <size_t I>				\
    decltype(auto) nestedExpr() ATTRIB		\
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
    
    /// Import assignement operator
    using Expr<T,TC,F>::operator=;
  };
  
  /// Combine the dynamic components of a tuple of dynamic comps, filling with each occurrence
  template <typename DcsOut,
	    typename..._DcsIn>
  CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
  auto dynamicCompsCombiner(const std::tuple<_DcsIn...>& dcsIns)
  {
    using DcsIns=
      std::tuple<_DcsIn...>;
    
    /// Result
    DcsOut dcsOut;
    
    UNROLL_FOR(int,i,0,2)
      {
	EXEC_FOR_ALL_TUPLE_IDS(IDcsIn,DcsIns,
			       
			       /// Input component on which we loop
			       using DcsIn=
			       std::tuple_element_t<IDcsIn,DcsIns>;
			       
			       /// List of dynamic components in common with result
			       using DcsCommonToOut=
			       TupleCommonTypes<DcsOut,DcsIn>;
			       
			       /// Value of all dynamic components
			       decltype(auto) dcsIn=
			       std::get<IDcsIn>(dcsIns);
				 
			       EXEC_FOR_ALL_TUPLE_IDS(IDcIn,DcsCommonToOut,
						      
						      const auto& dcIn=
						      std::get<IDcIn>(dcsIn);
						      
						      auto& dcOut=
						      std::get<IDcIn>(dcsOut);
						      
						      if(i==0)
							dcOut=dcIn;
						      else
							if(dcOut!=dcIn)
							  crash("unmatched dynamic comps among expressions");
						      ));
      }
    UNROLL_FOR_END;
    
    return
      dcsOut;
  }
}

#endif
