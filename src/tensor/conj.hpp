#ifndef _CONJ_HPP
#define _CONJ_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file conj.hpp

#include <metaProgramming/universalReference.hpp>
#include <tensor/unaryExpr.hpp>

namespace nissa
{
  DECLARE_COMPONENT(ComplId,int,2);
  
  /// Real component index
  constexpr inline ComplId Re=0;
  
  /// Imaginary component index
  constexpr inline ComplId Im=1;
  
  DEFINE_FEATURE(ConjugatorFeat);
  
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
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    EvalTo eval(const TD&...td) const
    {
      /// Compute the real or imaginary component
      const ComplId& reIm=
	std::get<ComplId>(std::make_tuple(td...));
      
      /// Nested result
      decltype(auto) nestedRes=
	this->nestedExpression.eval(td...);
      
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
      UnEx(FORWARD_MEMBER_VAR(Conjugator,oth,nestedExpression))
    {
    }
  };
  
  /// Takes the conjugate of e
  template <typename _E>
  auto conj(_E&& e,
	    UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR)
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
      
    return
      Conjugator<decltype(e),Comps,EvalTo>(std::forward<_E>(e));
  }
}

#endif
