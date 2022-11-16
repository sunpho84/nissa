#ifndef _CONJ_HPP
#define _CONJ_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file conj.hpp

#include <metaProgramming/universalReference.hpp>
#include <metaProgramming/dispatchStrategy.hpp>
#include <tensor/unaryExpr.hpp>

namespace nissa
{
  DECLARE_COMPONENT(ComplId,int,2);
  
#define FOR_REIM_PARTS(NAME)		\
  FOR_ALL_COMPONENT_VALUES(ComplId,NAME)
  
  /// Real component index - we cannot rely on a constexpr inline as the compiler does not propagate it correctly
#define Re ComplId(0)
  
  /// Imaginary component index
#define Im ComplId(1)
  
  /// Returns the real part, subscribing the complex component to Re value
  template <typename _E>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
  decltype(auto) real(_E&& e)
  {
    return
      //e(ComplId(0));
      e(Re);
  }
  
  /// Returns the imaginary part, subscribing the complex component to Im value
  template <typename _E>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
  decltype(auto) imag(_E&& e)
  {
    return
      e(Im);
  }
  
  /////////////////////////////////////////////////////////////////
  
  DEFINE_FEATURE(Conjugator);
  
#define THIS					\
  Conjugator<_E,_Comps,_EvalTo>
  
#define UNEX					\
    UnaryExpr<THIS,				\
	      _E,				\
	      _Comps,				\
	      _EvalTo>
  
  /// Conjugator of an expression
  template <typename _E,
	    typename _Comps,
	    typename _EvalTo>
  struct Conjugator :
    ConjugatorFeat<THIS>,
    UNEX
  {
    /// Import unary expression
    using UnEx=
      UNEX;
    
#undef UNEX
#undef THIS
    
    /// Components
    using Comps=
      _Comps;
    
    /// Type returned when evaluating
    using EvalTo=
      _EvalTo;
    
    /// Expression flags
    static constexpr ExprFlags Flags=
      unsetEvalToRef<setStoreByRefTo<false,UnEx::NestedFlags>>;
    
    /// Evaluate
    template <typename...TD>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    EvalTo eval(const TD&...td) const
    {
      /// Compute the real or imaginary component
      const ComplId& reIm=
	std::get<ComplId>(std::make_tuple(td...));
      
      /// Nested result
      decltype(auto) nestedRes=
	this->nestedExpr.eval(td...);
      
      if(reIm==0)
	return nestedRes;
      else
	return -nestedRes;
    }
    
    /// Construct
    template <typename C>
    Conjugator(C&& conjExpression) :
      UnEx(std::forward<C>(conjExpression))
    {
    }
    
    /// Move constructor
    Conjugator(Conjugator&& oth) :
      UnEx(FORWARD_MEMBER_VAR(Conjugator,oth,nestedExpr))
    {
    }
  };
  
  namespace internal
  {
    /// Decide the strategy to take a conjugate
    struct _ConjStrategy
    {
      /// Possible strategies
      enum{NO_COMPL,IS_CONJUGATOR,DEFAULT};
      
      /// Get the strategy for expression E
      template <typename E>
      using GetForExpr=
	std::integral_constant<int,
	(not tupleHasType<typename E::Comps,ComplId>)?
	NO_COMPL:
	    ((isConjugator<E>)?
	      IS_CONJUGATOR:
	     DEFAULT)>*;
      
      DECLARE_DISPATCH_STRATEGY(NoCompl,NO_COMPL);
      
      DECLARE_DISPATCH_STRATEGY(IsConjugator,IS_CONJUGATOR);
      
      DECLARE_DISPATCH_STRATEGY(Default,DEFAULT);
    };
    
    /// No conjugate is present
    template <typename _E>
    decltype(auto) _conj(_E&& e,
			 _ConjStrategy::NoCompl)
    {
      return
	e;
    }
    
    /// Returns the original expression
    template <typename _E>
    decltype(auto) _conj(_E&& e,
			 _ConjStrategy::IsConjugator)
    {
      return
	e.nestedExpr;
    }
    
    /// Takes the conjugate of e
    template <typename _E>
    auto _conj(_E&& e,
	       _ConjStrategy::Default)
    {
      /// Base expression
      using E=
	std::decay_t<_E>;
      
      /// Fundamental type
      using EvalTo=
	typename E::EvalTo;
      
      /// Components
      using Comps=
	typename E::Comps;
      
      // /// Position of complex component
      // static constexpr size_t complPos=
      // 	firstOccurrenceOfTypeInTuple<ComplId,Comps>;
      
      return
	Conjugator<decltype(e),Comps,EvalTo>(std::forward<_E>(e));
    }
  };
  
  /// Dispatch the correct conjugation
  template <typename _E>
  decltype(auto) conj(_E&& e)
  {
    return
      internal::_conj(std::forward<_E>(e),internal::_ConjStrategy::GetForExpr<std::decay_t<_E>>());
  }
}

#endif
