#ifndef _TRANSP_HPP
#define _TRANSP_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/nodes/transp.hpp

//#include <typeinfo>

//#include <routines/ios.hpp>

#include <expr/comp.hpp>
#include <expr/comps.hpp>
#include <expr/node.hpp>
#include <metaprogramming/universalReference.hpp>

namespace nissa
{
  PROVIDE_FEATURE(Transposer);
  
  /// Transposer
  ///
  /// Forward declaration to capture the components
  template <DerivedFromNode _E,
	    typename _Comps,
	    typename _Fund>
  struct Transposer;
  
#define THIS					\
  Transposer<_E,CompsList<C...>,_Fund>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// Transposer
  ///
  template <DerivedFromNode _E,
	    DerivedFromComp...C,
	    typename _Fund>
  struct THIS :
    TransposerFeat,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    /// Components
    using Comps=
      CompsList<C...>;
    
    /// Fundamental tye
    using Fund=_Fund;
    
    /// Transposed expression
    NodeRefOrVal<_E> transpExpr;
    
    /// Type of the transposed expression
    using TranspExpr=std::decay_t<_E>;
    
    /// Executes according to subexpr
    static constexpr ExecSpace execSpace=
      TranspExpr::execSpace;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    decltype(auto) getDynamicSizes() const
    {
      return transpExpr.getDynamicSizes();
    }
    
    /// Returns whether can assign
    INLINE_FUNCTION
    bool canAssign()
    {
      return transpExpr.canAssign();
    }
    
    /// This is a lightweight object
    static constexpr bool storeByRef=false;
    
    /// Import assignment operator
    using Base::operator=;
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=
      TranspExpr::canAssignAtCompileTime;
    
    /////////////////////////////////////////////////////////////////
    
    //// Returns a transposer on a different expression
    template <typename T>
    INLINE_FUNCTION
    decltype(auto) recreateFromExprs(T&& t) const
    {
      return transp(std::forward<T>(t));
    }
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)					\
    /*! Returns a reference */					\
    INLINE_FUNCTION						\
    auto getRef() ATTRIB					\
    {								\
      return transp(transpExpr.getRef());			\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
    /// Type obtained reinterpreting the fund
    template <typename NFund>
    using ReinterpretFund=
      Transposer<SameRefAs<_E,typename std::decay_t<_E>::template ReinterpretFund<NFund>>,
		 CompsList<C...>,
		 NFund>;
    
    /////////////////////////////////////////////////////////////////
    
    /// Evaluate
    template <typename...TD>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    Fund eval(const TD&...td) const
    {
      return transpExpr(transp(td)...);
    }
    
    /// Move construct
    INLINE_FUNCTION constexpr
    Transposer(Transposer&& oth) = default;
    
    /// Copy construct
    INLINE_FUNCTION constexpr
    Transposer(const Transposer& oth) = default;
    
    /// Construct
    template <typename T>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    Transposer(T&& arg)
      requires(std::is_same_v<std::decay_t<T>,std::decay_t<_E>>)
      : transpExpr{std::forward<T>(arg)}
    {
    }
  };
  
  /// Transpose an expression
  template <DerivedFromNode _E>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
  decltype(auto) transp(_E&& e)
  {
#if 0
    master_printf("Now inside transp\n");
#endif
    
    /// Base passed type
    using E=
      std::decay_t<_E>;
    
    if constexpr(isTransposer<E>)
      return e.transpExpr;
    else
      {
	using ArgComps=
	  typename E::Comps;
	
	/// Components
	using Comps=
	  TranspMatrixTensorComps<ArgComps>;
	
	if constexpr(not compsAreTransposable<Comps>)
	  {
#if 0
	    master_printf("no need to transpose, returning the argument, which is %p %s %s",&e,typeid(_E).name(),(std::is_lvalue_reference_v<decltype(e)>?"&":(std::is_rvalue_reference_v<decltype(e)>?"&&":"")));
#endif
	    
	    return RemoveRValueReference<_E>(e);
	  }
	else
	  {
	    /// Type returned when evaluating the expression
	    using Fund=
	      typename E::Fund;
	    
	    return
	      Transposer<decltype(e),Comps,Fund>(std::forward<_E>(e));
	  }
      }
  }
}

#endif
