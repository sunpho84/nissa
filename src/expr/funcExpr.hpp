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
    template <typename...Ci>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    Fund eval(const Ci&...c) const
    {
      return func(c...);
    }
    
    /// Construct
    template <DerivedFromComp...TD>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
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
		execOnCPUAndGPU;
    
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    FuncNodeWrapper(F f) :
      f(f)
    {
    }
    
    /// Evaluate
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    auto operator()(const C&...c) const
    {
      return f(c...);
    }
  };
  
  template <typename C,
	    ExecSpace ES=execOnCPUAndGPU,
	    typename F,
	    DerivedFromComp...TD>
  FuncExpr<FuncNodeWrapper<F,C,ES>,C,decltype(std::apply(std::declval<F>(),std::declval<C>()))> funcNodeWrapper(const F& f,
														const CompsList<TD...>& ds)
  {
    return {FuncNodeWrapper<F,C,ES>(f),ds};
  }
}

#endif
