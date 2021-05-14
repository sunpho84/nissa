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
  DECLARE_COMPONENT(Complex,int,2);
  
  /// Real component index
  constexpr inline Complex Re=0;
  
  /// Imaginary component index
  constexpr inline Complex Im=1;
  
  DEFINE_FEATURE(ConjFeat);
  
#define THIS					\
  Conj<_E,_Comps,_Fund>
  
#define UNEX					\
    UnaryExpr<THIS,				\
	      _E,				\
	      _Comps,				\
	      _Fund>
  
  /// Conjugator of an expression
  template <typename _E,
	    typename _Comps,
	    typename _Fund>
  struct Conj :
    ConjFeat<THIS>,
    UNEX
  {
    using UnEx=
      UNEX;
    
#undef UNEX
#undef THIS
    
    /// Components
    using Comps=
      _Comps;
    
    /// Fundamental type of the expression
    using Fund=
      _Fund;
    
    /// Expression flags
    static constexpr ExprFlags Flags=
      unsetEvalToRef<setStoreByRefTo<false,UnEx::NestedFlags>>;
    
    /// Evaluate
    template <typename...TD>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    Fund eval(const TD&...td) const
    {
      const Complex& reIm=
	std::get<Complex>(std::make_tuple(td...));
      
      decltype(auto) val=
	this->nestedExpression.eval(td...);
      
      if(reIm==0)
	return val;
      else
	return -val;
    }
    
    /// Construct
    template <typename C>
    Conj(C&& conjExpression) :
      UnEx(std::forward<C>(conjExpression))
    {
    }
    
    /// Move constructor
    Conj(Conj&& oth) :
      UnEx(FORWARD_MEMBER_VAR(Conj,oth,nestedExpression))
    {
    }
  };
  
  /// Takes the conjugate of e
  template <typename _E>
  auto conj(_E&& e,
	    UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR)
  {
    using E=
      std::decay_t<_E>;
    
    using Fund=
      typename E::Fund;
    
    using Comps=
      typename E::Comps;
      
    return
      Conj<decltype(e),Comps,Fund>(std::forward<_E>(e));
  }
}

#endif
