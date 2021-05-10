#ifndef _CONJ_HPP
#define _CONJ_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file conj.hpp

#include <metaProgramming/universalReference.hpp>
#include <tensor/expr.hpp>
#include <tensor/refCatcher.hpp>

namespace nissa
{
  DECLARE_COMPONENT(Complex,int,2);
  
  /// Real component index
  constexpr inline Complex Re=0;
  
  /// Imaginary component index
  constexpr inline Complex Im=1;
  
  /// Conjugator of an expression
  template <typename T,
	    ExprFlags _Flags>
  struct Conj : Expr<Conj<T,_Flags>,
		     typename T::Comps,
		     typename T::Fund,
		     unsetEvalToRef<setStoreByRefTo<false,_Flags>>>
  {
    /// Components
    using Comps=
      typename T::Comps;
    
    /// Fundamental type of the expression
    using Fund=
      typename T::Fund;
    
    /// Type of the bound expression
    using BoundExpression=
      ConditionalRef<getStoreByRef<_Flags>
      ,ConditionalConst<getEvalToConst<_Flags>,T>>;
    
    /// Bound expression
    BoundExpression boundExpression;
    
    /// Evaluate
    template <typename...TD>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    Fund eval(const TD&...td) const
    {
      const Complex& reIm=
	std::get<Complex>(std::make_tuple(td...));
      
      decltype(auto) val=boundExpression.eval(td...);
      
      if(reIm==0)
	return val;
      else
	return -val;
    }
    
    /// Construct
    template <typename C>
    Conj(C&& boundExpression) :
      boundExpression(std::forward<C>(boundExpression))
    {
    }
    
    /// Move constructor
    Conj(Conj&& oth) :
      boundExpression(FORWARD_MEMBER_VAR(Conj,oth,boundExpression))
    {
    }
  };
  
  template <typename _E>
  auto conj(_E&& e,
	    UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR)
  {
    using CH=
      RefCatcherHelper<_E,decltype(e)>;
    
    return Conj<typename CH::E,CH::Flags>(std::forward<_E>(e));
  }
  
}

#endif
