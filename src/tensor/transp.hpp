#ifndef _TRANSP_HPP
#define _TRANSP_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file transp.hpp

#include <metaProgramming/universalReference.hpp>
#include <tensor/unaryExpr.hpp>

namespace nissa
{
  
#define THIS					\
  Transp<_E,_Comps,_Fund>

#define UNEX					\
    UnaryExpr<THIS,				\
	      _E,				\
	      _Comps,				\
	      _Fund>
  
  /// Transposer of an expression
  template <typename _E,
	    typename _Comps,
	    typename _Fund>
  struct Transp : UNEX
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
      unsetStoreByRef<UnEx::NestedFlags>;
    
#define DECLARE_EVAL(ATTRIB)				\
    /*! Evaluate */					\
    template <typename...TD>				\
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr		\
    decltype(auto) eval(const TD&...td) ATTRIB		\
    {							\
      return						\
	this->nestedExpression(td.transp()...);		\
    }
    
    DECLARE_EVAL(const);
    
    DECLARE_EVAL(/*non const*/);
    
#undef DECLARE_EVAL
    
    /// Construct
    template <typename C>
    Transp(C&& boundExpression) :
      UnEx(std::forward<C>(boundExpression))
    {
    }
    
    /// Move constructor
    Transp(Transp&& oth) :
      UnEx(FORWARD_MEMBER_VAR(Transp,oth,boundExpression))
    {
    }
  };
  
  template <typename _E>
  auto transp(_E&& e,
	      UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR)
  {
    using E=
      std::decay_t<_E>;
    
    using Fund=
      typename E::Fund;
    
    using Comps=
      TransposeTensorComps<typename E::Comps>;
    
    return
      Transp<decltype(e),Comps,Fund>(std::forward<_E>(e));
  }
}

#endif
