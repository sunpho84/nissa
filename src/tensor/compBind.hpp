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

namespace nissa
{
  template <typename E,
	    typename TC,
	    typename BC,
	    typename IC,
	    ExprFlags _Flags>
  struct _CompBinderFactory;
  
  template <typename E,
	    typename...TCs,
	    typename...BCs,
	    size_t...ICs,
	    ExprFlags _Flags>
  struct _CompBinderFactory<E,
			    TensorComps<TCs...>,
			    TensorComps<BCs...>,
			    std::index_sequence<ICs...>,
			    _Flags>
  {
    using BoundExpressionComponents=
      TensorComps<TCs...>;
    
    using _BoundExpression=
      ConditionalRef<getStoreByRef<_Flags>
      ,ConditionalConst<getEvalToConst<_Flags>,E>>;
    
    using _BoundComponents=
      TensorComps<BCs...>;
    
    using BoundFund=
      typename E::Fund;
    
    using FilteredComponents=
      TupleFilterAllTypes<BoundExpressionComponents,
			  _BoundComponents>;
    
    struct _CompBinder :
      Expr<_CompBinder,FilteredComponents,BoundFund,setStoreByRefTo<false,_Flags>>
    {
      using Comps=
	FilteredComponents;
      
      using BoundComponents=
	_BoundComponents;
      
      using BoundExpression=
	_BoundExpression;
      
      /// Fundamental type of the expression
      using Fund=
	BoundFund;
      
      /// Fundamental type, possibly with constant attribute if required via flags
      using MaybeConstFund=
	ConditionalConst<getEvalToConst<_Flags>,const Fund>;
      
      /// Possibly constant fundamental type, possibly referencing
      using MaybeConstMaybeRefToFund=
	ConditionalRef<getEvalToRef<_Flags>,MaybeConstFund>;
      
      /// Bound expression
      BoundExpression boundExpression;
      
      /// Components that have been bound
      BoundComponents boundComponents;
      
      /// Constant access, returning constant reference or type
      template <typename...TD>
      CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
      const MaybeConstMaybeRefToFund eval(const TD&...td) const
      {
	return boundExpression.eval(std::get<ICs>(boundComponents)...,td...);
      }
      
      /// Non const access, possibly still returning a constant reference
      template <typename...TD>
      CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
      MaybeConstMaybeRefToFund eval(const TD&...td)
      {
	return boundExpression.eval(std::get<ICs>(boundComponents)...,td...);
      }
      
      /// Construct
      template <typename T>
      _CompBinder(T&& boundExpression,
		  const BoundComponents& boundComponents) :
	boundExpression(std::forward<T>(boundExpression)),
	boundComponents(boundComponents)
      {
      }
      
      /// Move constructor
      _CompBinder(_CompBinder&& oth) :
	boundExpression(FORWARD_MEMBER_VAR(_CompBinder,oth,boundExpression)),
	boundComponents(oth.boundComponents)
      {
      }
      
      /// Copy constructor
      _CompBinder(const _CompBinder& oth) :
	boundExpression(oth.boundExpression),
	boundComponents(oth.boundComponents)
      {
	static_assert(std::is_lvalue_reference_v<BoundExpression>,"Makes no sense");
      }
    };
  };
  
#define _CBF					\
  _CompBinderFactory<E,				\
		     typename E::Comps,		\
		     BC,						\
		     std::make_index_sequence<std::tuple_size_v<BC>>,	\
		     Flags>::_CompBinder
  
  template <typename E,
	    typename BC,
	    ExprFlags Flags>
  struct CompBinder : _CBF
  {
    using CBF=
      typename _CBF;
    
#undef _CBF
    
    using CBF::CBF;
    
    using BoundComponents=
      typename CBF::BoundComponents;
    
    using BoundExpression=
      typename CBF::BoundExpression;
    
    using Fund=
      typename CBF::Fund;
  };
  
  template <typename _E,
	    typename ACTUAL_TYPE> /// corresponding to decltype(e)
  struct _CompBinderTraits
  {
    using E=
      std::remove_const_t<std::decay_t<_E>>;
    
    static constexpr bool storeByRef=
		std::is_lvalue_reference_v<ACTUAL_TYPE>;
    
    static constexpr bool needsToBeMoveConstructed=
		std::is_rvalue_reference_v<ACTUAL_TYPE>;
    
    static constexpr bool isVal=
		not std::is_reference_v<ACTUAL_TYPE>;
    
    static constexpr bool canBeMoveConstructed=
		std::is_move_constructible_v<E>;
    
    static constexpr bool canBeCopyConstructed=
		std::is_copy_constructible_v<E>;
    
    static constexpr bool passAsConst=
		std::is_const_v<std::remove_reference_t<ACTUAL_TYPE>>;
    
    static_assert(canBeMoveConstructed or not needsToBeMoveConstructed,"Would need to move-construct, but the move constructor is not available");
    
    static_assert(canBeCopyConstructed or not isVal,
		  "Would need to copy-construct, but the copy constructor is not available or the inner object must be stored by ref");
    
    static_assert((not getStoreByRef<E::Flags>) or not isVal,
		  "Would need to store by val, but the inner object flags indicated to store by ref");
    
    static constexpr ExprFlags Flags
		=setStoreByRefTo<storeByRef
		,addEvalToConstIf<passAsConst,E::Flags>>;
  };
  
  template <typename _E,
	    typename...BCs>
  auto compBind(_E&& e,
		const TensorComps<BCs...>& bc,
		UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR)
  {
    using CH=_CompBinderTraits<_E,decltype(e)>;
    
    return
      CompBinder<typename CH::E,TensorComps<BCs...>,CH::Flags>(std::forward<_E>(e),bc);
  }
  
  template <typename E,
	    typename BC,
	    ExprFlags Flags,
	    typename...BCs>
  auto compBind(const CompBinder<E,BC,Flags>& cb,
		const TensorComps<BCs...>& bcs)
  {
    return compBind(cb.boundExpression,std::tuple_cat(cb.boundComponents,bcs));
  }
}

#endif
