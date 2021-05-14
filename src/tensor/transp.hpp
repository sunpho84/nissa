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
  DEFINE_FEATURE(TransposerFeat);
  
#define THIS					\
  Transposer<_E,_Comps,_EvalTo>

#define UNEX					\
    UnaryExpr<THIS,				\
	      _E,				\
	      _Comps,				\
	      _EvalTo>
  
  /// Transposeroser of an expression
  template <typename _E,
	    typename _Comps,
	    typename _EvalTo>
  struct Transposer :
    TransposerFeat<THIS>,
    UNEX
  {
    /// Import the unary expression
    using UnEx=
      UNEX;
    
#undef UNEX
#undef THIS
    
    /// Components
    using Comps=
      _Comps;
    
    /// Fundamental type of the expression
    using EvalTo=
      _EvalTo;
    
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
    Transposer(C&& boundExpression) :
      UnEx(std::forward<C>(boundExpression))
    {
    }
    
    /// Move constructor
    Transposer(Transposer&& oth) :
      UnEx(FORWARD_MEMBER_VAR(Transposer,oth,boundExpression))
    {
    }
  };
  
  /// Returns the transposer of e
  template <typename _E>
  auto transp(_E&& e,
	      UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR)
  {
    /// Decayed type
    using E=
      std::decay_t<_E>;
    
    /// Type returned when evaluating
    using EvalTo=
      typename E::EvalTo;
    
    /// Components
    using Comps=
      TransposeTensorComps<typename E::Comps>;
    
    return
      Transposer<decltype(e),Comps,EvalTo>(std::forward<_E>(e));
  }
}

#endif
