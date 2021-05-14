#ifndef _COMPBIND_HPP
#define _COMPBIND_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file compBind.hpp

#include <metaProgramming/refOrVal.hpp>
#include <metaProgramming/universalReference.hpp>
#include <tensor/component.hpp>
#include <tensor/expr.hpp>
#include <tensor/unaryExpr.hpp>

namespace nissa
{
  DEFINE_FEATURE(CompBindFeat);
  
  /// Component binder
  ///
  /// Forward declaration to capture the index sequence
  template <typename _BC,
	    typename _IC,
	    typename _E,
	    typename _Comps,
	    typename _Fund>
  struct CompBinder;
  
#define THIS					\
  CompBinder<_BC,std::index_sequence<ICs...>,_E,_Comps,_Fund>

#define UNEX					\
    UnaryExpr<THIS,				\
	      _E,				\
	      _Comps,				\
	      _Fund>
  
  /// Component binder
  ///
  template <typename _BC,
	    size_t...ICs,
	    typename _E,
	    typename _Comps,
	    typename _Fund>
  struct THIS :
    CompBindFeat<THIS>,
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
    using Fund=
      _Fund;
    
    /// Expression flags
    static constexpr ExprFlags Flags=
      UnEx::NestedFlags;
    
    /// List of bound components
    using BoundComps=
      _BC;
    
    /// Components that have been bound
    BoundComps boundComps;
    
#define DECLARE_EVAL(ATTRIB)						\
									\
    /*! Constant access, returning ATTRIB reference or type */		\
    template <typename...TD>						\
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr				\
    decltype(auto) eval(const TD&...td) ATTRIB				\
    {									\
      return								\
	this->nestedExpression.eval(std::get<ICs>(boundComps)...,td...); \
    }
    
    DECLARE_EVAL(const);
    
    DECLARE_EVAL(/*non const*/);
    
#undef DECLARE_EVAL
    
    /// Construct
    template <typename T>
    CompBinder(T&& boundExpression,
	       const BoundComps& boundComps) :
      UnEx(std::forward<T>(boundExpression)),
      boundComps(boundComps)
    {
    }
    
    /// Move constructor
    CompBinder(CompBinder&& oth) :
      UnEx(FORWARD_MEMBER_VAR(CompBinder,oth,nestedExpression)),
      boundComps(oth.boundComps)
    {
    }
    
    /// Copy constructor
    CompBinder(const CompBinder& oth) :
      UnEx(oth.nestedExpression),
      boundComps(oth.boundComps)
    {
      static_assert(std::is_lvalue_reference_v<typename UnEx::BoundExpression>,"It makes no sense");
    }
  };
  
  template <typename _E,
	    typename...BCs>
  auto compBind(_E&& e,
		const TensorComps<BCs...>& bc,
		UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR)
  {
    /// Base passed type
    using E=
      std::decay_t<_E>;
    
    /// Fundamental type of the expression
    using Fund=
      typename E::Fund;
    
    /// Components to bind
    using BoundComps=
      TensorComps<BCs...>;
    
    /// Visible components
    using Comps=
      TupleFilterAllTypes<typename E::Comps,
			  BoundComps>;
    return
      CompBinder<BoundComps,
		 std::make_index_sequence<sizeof...(BCs)>,
		 decltype(e),
		 Comps,
		 Fund>(std::forward<_E>(e),bc);
  }
  
  template <typename CB,
	    typename...BCs>
  auto compBind(const CompBindFeat<CB>& cb,
		const TensorComps<BCs...>& bcs)
  {
    return
      compBind(cb.defeat().nestedExpression,
	       std::tuple_cat(cb.deFeat().boundComps,bcs));
  }
}

#endif
