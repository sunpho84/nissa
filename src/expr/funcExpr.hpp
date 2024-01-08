#ifndef _FUNCEXPR_HPP
#define _FUNCEXPR_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/nodes/funcExpr.hpp

#include <expr/comps.hpp>
#include <expr/dynamicCompsProvider.hpp>
#include <expr/execSpace.hpp>
#include <expr/nodeDeclaration.hpp>
#include <expr/subExprs.hpp>
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
    FuncExprFeat,
    NoSubExprs,
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
    
    /// Type of the function
    using Func=_Func;
    
    /// Executes where function can
    static constexpr ExecSpace execSpace=
      Func::execSpace;
    
    /// Function to be called
    const Func func;
    
    /// Sizes of the dynamic components
    DynamicComps dynamicSizes;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
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
    
    /// Type obtained reinterpreting the fund
    template <typename NFund>
    using ReinterpretFund=
      FuncExpr<_Func,CompsList<C...>,NFund>;
    
    /// Evaluate
    template <typename...Ci>
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    Fund eval(const Ci&...c) const
    {
      return func(c...);
    }
    
    /// Construct
    template <DerivedFromComp...TD>
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    FuncExpr(const Func& func,
	     const CompsList<TD...>& td) :
      func(func),
      dynamicSizes(tupleGetSubset<DynamicComps>(td))
    {
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  /// Wraps a function callable with comps C
  template <typename F,
	    typename C,
	    ExecSpace ES=execOnCPUAndGPU>
  struct FuncNodeWrapper;
  
  /// Wraps a function callable with comps C
  template <typename F,
	    DerivedFromComp...C,
	    ExecSpace ES>
  requires(std::is_invocable_v<F,C...>)
  struct FuncNodeWrapper<F,CompsList<C...>,ES>
  {
    F f;
    
    /// Can run on both GPU and CPU as it is trivially copyable
    static constexpr ExecSpace execSpace=
      ES;
    
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    FuncNodeWrapper(F f) :
      f(f)
    {
    }
    
    /// Evaluate
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    auto operator()(const C&...c) const
    {
      return f(c...);
    }
  };
  
  /// Creates a wrapper around a callable
  template <typename C,
	    ExecSpace ES=ExecSpace{execOnCPUAndGPU},
	    typename F,
	    typename Fund=decltype(std::apply(std::declval<F>(),std::declval<C>())),
	    DerivedFromComp...TD>
  constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
  auto funcNodeWrapper(F&& f,
		       const CompsList<TD...>& ds)
  {
    return FuncExpr<FuncNodeWrapper<F,C,ES>,C,Fund>(FuncNodeWrapper<F,C,ES>(std::forward<F>(f)),ds);
  }
}

#endif
