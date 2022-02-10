#ifndef _TRACE_HPP
#define _TRACE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file trace.hpp

#include <tensor/loopOnAllComponents.hpp>
#include <metaProgramming/universalReference.hpp>
#include <tensor/unaryExpr.hpp>

namespace nissa
{
  namespace internal
  {
    /// Classifies the components, determining which one are visible or traced
    ///
    /// Internal implementation, forward declaration
    template <typename TC>
    struct _TraceCompsComputer;
    
    /// Classifies the components, determining which one are visible or traced
    template <typename...Tc>
    struct _TraceCompsComputer<TensorComps<Tc...>>
    {
      template <typename T>
      struct Classify
      {
	static constexpr RwCl rC=
	   T::rC;
	
	static constexpr bool isTraced=
	  T::isMatrix and
	  (std::is_same_v<T,Tc>||...) and
	  (std::is_same_v<typename T::Transp,Tc>||...);
	
	///
	using VisiblePart=
	  std::conditional_t<isTraced,TensorComps<>,TensorComps<T>>;
	
	using TracedPart=
	  std::conditional_t<isTraced and rC==ROW,TensorComps<T>,TensorComps<>>;
      };
      
      using VisibleComps=
	TupleCat<typename Classify<Tc>::VisiblePart...>;
      
      using TracedComps=
	TupleCat<typename Classify<Tc>::TracedPart...>;
    };
  }
  
  DEFINE_FEATURE(Tracer);
  
#define THIS					\
  Tracer<_TracedComps,_E,_Comps,_EvalTo>
  
#define UNEX					\
    UnaryExpr<THIS,				\
	      _E,				\
	      _Comps,				\
	      _EvalTo>
  
  /// Transposeroser of an expression
  template <typename _TracedComps,
	    typename _E,
	    typename _Comps,
	    typename _EvalTo>
  struct Tracer :
    TracerFeat<THIS>,
    UNEX
  {
    /// Import the unary expression
    using UnEx=
      UNEX;
    
#undef UNEX
#undef THIS
    
    /// Traced comps
    using TracedComps=
      _TracedComps;
    
    /// Components
    using Comps=
      _Comps;
    
    /// Fundamental type of the expression
    using EvalTo=
      _EvalTo;
    
    /// Expression flags
    static constexpr ExprFlags Flags=
      unsetStoreByRef<UnEx::NestedFlags>;
    
    ///Evaluate
    template <typename...NTCs>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    EvalTo eval(const NTCs&...nTCs) const
    {
      /// Result
      EvalTo res(0);
      
      loopOnAllComponents<TracedComps>(TensorComps<>(),
				       [this,&res,&nTCs...](const auto&...tCs) INLINE_ATTRIBUTE
				       {
					 /// First argument
					 res+=
					   this->nestedExpr(nTCs...,tCs...,tCs.transp()...);
				       });
      
      return
	res;
    }
    
    /// Construct
    template <typename C>
    Tracer(C&& boundExpression) :
      UnEx(std::forward<C>(boundExpression))
    {
    }
    
    /// Move constructor
    Tracer(Tracer&& oth) :
      UnEx(FORWARD_MEMBER_VAR(Tracer,oth,nestedExpr))
    {
    }
  };
  
  /// Returns the tracer of e
  template <typename _E,
	    UNPRIORITIZE_DEFAULT_VERSION_TEMPLATE_PARS>
  auto trace(_E&& e,
	      UNPRIORITIZE_DEFAULT_VERSION_ARGS)
  {
    UNPRIORITIZE_DEFAULT_VERSION_ARGS_CHECK;
    
    /// Decayed type
    using E=
      std::decay_t<_E>;
    
    /// Type returned when evaluating
    using EvalTo=
      typename E::EvalTo;
    
    /// Computes the trace components
    using TCC=
      internal::_TraceCompsComputer<typename E::Comps>;
    
    /// Components
    using Comps=
      typename TCC::VisibleComps;
    
    /// Traced Components
    using TracedComps=
      typename TCC::TracedComps;
    
    return
      Tracer<TracedComps,decltype(e),Comps,EvalTo>(std::forward<_E>(e));
  }
}

#endif
