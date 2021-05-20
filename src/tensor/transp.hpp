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
  DEFINE_FEATURE(Transposer);
  
#define THIS					\
  Transposer<_E,_Comps,_CompsMeldBarriers,_EvalTo>

#define UNEX					\
    UnaryExpr<THIS,				\
	      _E,				\
	      _Comps,				\
	      _CompsMeldBarriers,		\
	      _EvalTo>
  
  /// Transposeroser of an expression
  template <typename _E,
	    typename _Comps,
	    typename _CompsMeldBarriers,
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
    
    /// Barrier to meld components
    using CompsMeldBarriers=
      _CompsMeldBarriers;
    
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
	this->nestedExpr(td.transp()...);		\
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
      UnEx(FORWARD_MEMBER_VAR(Transposer,oth,nestedExpr))
    {
    }
  };
  
  /// Returns the transposer of e
  template <typename _E,
	    UNPRIORITIZE_DEFAULT_VERSION_TEMPLATE_PARS>
  auto transp(_E&& e,
	      UNPRIORITIZE_DEFAULT_VERSION_ARGS)
  {
    UNPRIORITIZE_DEFAULT_VERSION_ARGS_CHECK;
    
    /// Decayed type
    using E=
      std::decay_t<_E>;
    
    /// Type returned when evaluating
    using EvalTo=
      typename E::EvalTo;
    
    /// Components
    using Comps=
      TransposeTensorComps<typename E::Comps>;
    
    /// Barriers between fusable components, obtained faking transposition of comps
    using CompsMeldBarriers=
      GetTensorCompsMeldBarriersFromPureRemapping<TransposeTensorComps<typename E::Comps>,typename E::Comps,typename E::CompsMeldBarriers>;
    
    return
      Transposer<decltype(e),Comps,CompsMeldBarriers,EvalTo>(std::forward<_E>(e));
  }
  
  /// Returns the transposer of e
  template <typename _E,
	    ENABLE_THIS_TEMPLATE_IF(isTransposer<_E>)>
  auto transp(_E&& e)
  {
    return
      e.nestedExpr;
  }
}

#endif
