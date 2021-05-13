#ifndef _TRANSP_HPP
#define _TRANSP_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file transp.hpp

#include <metaProgramming/universalReference.hpp>
#include <tensor/expr.hpp>
#include <tensor/refCatcher.hpp>

namespace nissa
{
  /// Transposer of an expression
  template <typename T,
	    ExprFlags _Flags>
  struct Transp : Expr<Transp<T,_Flags>,
		       TransposeTensorComps<typename T::Comps>,
		       typename T::Fund,
		       unsetEvalToRef<unsetStoreByRef<_Flags>>>
  {
    /// Components
    using Comps=
      TransposeTensorComps<typename T::Comps>;
    
    /// Fundamental type of the expression
    using Fund=
      typename T::Fund;
    
    /// Type of the bound expression
    using BoundExpression=
      ConditionalRef<getStoreByRef<_Flags>
      ,ConditionalConst<getEvalToConst<_Flags>,T>>;
    
    /// Bound expression
    BoundExpression boundExpression;
    
    /// Dynamic sizes
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) getDynamicSizes() const
    {
      return boundExpression.getDynamicSizes();
    }
    
    /// Evaluate
    template <typename...TD>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    Fund eval(const TD&...td) const
    {
      return boundExpression(td.transp()...);
    }
    
    /// Construct
    template <typename C>
    Transp(C&& boundExpression) :
      boundExpression(std::forward<C>(boundExpression))
    {
    }
    
    /// Move constructor
    Transp(Transp&& oth) :
      boundExpression(FORWARD_MEMBER_VAR(Transp,oth,boundExpression))
    {
    }
  };
  
  template <typename _E>
  auto transp(_E&& e,
	      UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR)
  {
    using CH=
      RefCatcherHelper<_E,decltype(e)>;
    
    return Transp<typename CH::E,CH::Flags>(std::forward<_E>(e));
  }
}

#endif
