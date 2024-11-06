#ifndef _TRACE_HPP
#define _TRACE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/trace.hpp

#include <expr/comps.hpp>
#include <expr/dynamicCompsProvider.hpp>
#include <expr/node.hpp>
#include <metaprogramming/arithmeticOperatorsViaCast.hpp>
#include <tuples/tupleCat.hpp>

namespace nissa
{
  /// Classifies the components, determining which one are visible or traced
  ///
  /// Internal implementation, forward declaration
  template <typename TC>
  struct TracerCompsDeducer;
  
  /// Classifies the components, determining which one are visible or traced
  template <typename...Tc>
  struct TracerCompsDeducer<CompsList<Tc...>>
  {
    template <typename T>
    struct Classify
    {
      static constexpr bool isTraced=
	isTransposableComp<T> and
	(std::is_same_v<T,Tc> or...) and
	(std::is_same_v<typename T::Transp,Tc> or...);
      
      ///
      using VisiblePart=
	std::conditional_t<isTraced,std::tuple<>,std::tuple<T>>;
      
      using TracedPart=
	std::conditional_t<isTraced and isRow<T>,std::tuple<T>,std::tuple<>>;
    };
    
    using VisibleComps=
      TupleCat<typename Classify<Tc>::VisiblePart...>;
    
    using TracedComps=
      TupleCat<typename Classify<Tc>::TracedPart...>;
  };
  
  /////////////////////////////////////////////////////////////////
  
  PROVIDE_FEATURE(Tracer);
  
  /// Tracer
  ///
  /// Forward declaration to capture the components
  template <typename Tc,
	    DerivedFromNode _E,
	    typename _Comps,
	    typename _Fund>
  struct Tracer;
  
#define THIS					\
  Tracer<CompsList<Tc...>,_E,CompsList<C...>,_Fund>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// Tracer
  ///
  template <DerivedFromComp...Tc,
	    DerivedFromNode _E,
	    DerivedFromComp...C,
	    typename _Fund>
  struct THIS :
    TracerFeat,
    SingleSubExpr<THIS>,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    /// Components
    using TracedComps=
      CompsList<Tc...>;
    
    /// Components
    using Comps=
      CompsList<C...>;
    
    /// Fundamental tye
    using Fund=_Fund;
    
    /// Expression to trace
    NodeRefOrVal<_E> subExpr;
    
    /// Type of the traced expression
    using TracedExpr=std::decay_t<_E>;
    
    /// Executes where traced reference
    static constexpr ExecSpace execSpace=
      TracedExpr::execSpace;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    decltype(auto) getDynamicSizes() const
    {
      return subExpr.getDynamicSizes();
    }
    
    /// Returns whether can assign
    INLINE_FUNCTION
    bool canAssign()
    {
      return false;
    }
    
    /// This is a lightweight object
    static constexpr bool storeByRef=false;
    
    /// Import assignment operator
    using Base::operator=;
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=false;
    
    /////////////////////////////////////////////////////////////////
    
    //// Returns a tracer on a different expression
    template <typename T>
    INLINE_FUNCTION
    decltype(auto) recreateFromExprs(T&& t) const
    {
      return trace(std::forward<T>(t));
    }
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)					\
    /*! Returns a reference */					\
    INLINE_FUNCTION						\
    auto getRef() ATTRIB					\
    {								\
      return trace(subExpr.getRef());			\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
    /// Type obtained reinterpreting the fund
    template <typename NFund>
    using ReinterpretFund=
      Tracer<CompsList<Tc...>,
	     SameRefAs<_E,typename std::decay_t<_E>::template ReinterpretFund<NFund>>,
	     CompsList<C...>,
	     NFund>;
    
    /////////////////////////////////////////////////////////////////
    
    /// Evaluate
    template <typename...NTc>
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    Fund eval(const NTc&...nTCs) const
    {
      if constexpr(std::tuple_size_v<TracedComps>)
	{
	  /// Result
	  Fund res;
	  setToZero(res);
	  
	  compsLoop<TracedComps>([this,&res,&nTCs...](const auto&...tCs) INLINE_ATTRIBUTE
	  {
	    /// First argument
	    res+=
	      this->subExpr(nTCs...,tCs...,transp(tCs)...);
	  },getDynamicSizes());
	  
	  return
	    res;
	}
      else
	return this->subExpr(nTCs...);
    }
    
    /// Construct
    template <typename T>
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    Tracer(T&& arg)
      requires(std::is_same_v<std::decay_t<T>,std::decay_t<_E>>)
      : subExpr(std::forward<T>(arg))
    {
    }
    
    // /// Copy constructor
    // HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    // Tracer(const Tracer& oth) :
    //   subExpr(oth.subExpr)
    // {
    // }
  };
  
  /// Trace an expression
  template <DerivedFromNode _E>
  HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
  decltype(auto) trace(_E&& e)
  {
    /// Base passed type
    using E=
      std::decay_t<_E>;
    
    using CompsDeducer=
      TracerCompsDeducer<typename E::Comps>;
    
    using TracedComps=typename CompsDeducer::TracedComps;
    
    using Comps=typename CompsDeducer::VisibleComps;
    
    using Fund=typename E::Fund;
    
    return
      Tracer<TracedComps,decltype(e),Comps,Fund>(std::forward<_E>(e));
  }
  
  /// Trace an expression
  template <DerivedFromTransposableComp...TracedComp,
	    DerivedFromNode _E>
  HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
  decltype(auto) traceOver(_E&& e)
  {
    /// Base passed type
    using E=
      std::decay_t<_E>;
    
    using TracedComps=
      std::tuple<TracedComp...>;
    
    using Comps=
      TupleFilterAllTypes<typename E::Comps,std::tuple<TracedComp...,typename TracedComp::Transp...>>;
    
    using Fund=typename E::Fund;
    
    return
      Tracer<TracedComps,decltype(e),Comps,Fund>(std::forward<_E>(e));
  }
}

#endif
