#ifndef _FUNCEXPR_HPP
#define _FUNCEXPR_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/nodes/funcExpr.hpp

#include <expr/comps.hpp>
#include <expr/dynamicCompsProvider.hpp>
#include <expr/nodeDeclaration.hpp>
#include <metaprogramming/universalReference.hpp>

namespace nissa
{
  PROVIDE_FEATURE(FuncExpr);
  
  /// Functional expression
  ///
  /// Forward declaration to capture the components
  template <typename Func,
	    typename _Comps,
	    typename _Fund>
  struct FuncExpr;
  
#define THIS					\
  FuncExpr<_Func,CompsList<C...>,_Fund>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// Functional expression
  ///
  /// Forward declaration
  template <typename _Func,
	    typename...C,
	    typename _Fund>
  struct FuncExpr<_Func,CompsList<C...>,_Fund> :
    FuncExprFeat<THIS>,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    /// List of dynamic comps
    using DynamicComps=
      typename DynamicCompsProvider<CompsList<C...>>::DynamicComps;
    
    /// Components
    using Comps=
      CompsList<C...>;
    
    /// Fundamental tye
    using Fund=_Fund;
    
    // /// Executes where allocated
    // static constexpr ExecSpace execSpace=
    //   SubNode<0>::execSpace;
    
    /// Type of the function
    using Func=_Func;
    
    /// Function to be called
    const Func func;
    
    /// Sizes of the dynamic components
    DynamicComps dynamicSizes;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    decltype(auto) getDynamicSizes() const
    {
      return dynamicSizes;
    }
    
    /// Returns whether can assign
    INLINE_FUNCTION
    constexpr bool canAssign()
    {
      return false;
    }
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=false;
    
    /// This is a lightweight object, hopefully
    static constexpr bool storeByRef=false;
    
    /////////////////////////////////////////////////////////////////
    
    // //// Returns a conjugator on a different expression
    // template <typename T>
    // INLINE_FUNCTION
    // decltype(auto) recreateFromExprs(T&& t) const
    // {
    //   return conj(std::forward<T>(t));
    // }
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)					\
    /*! Returns a reference */					\
    INLINE_FUNCTION						\
    auto getRef() ATTRIB					\
    {								\
      return *this;						\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /// Evaluate
    template <typename...TD>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    Fund eval(const TD&...td) const
    {
      return func(td...);
    }
    
    /// Construct
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    FuncExpr(const Func& func) :
      func(func)
    {
    }
  };
}

#endif
