#ifndef _COMPBIND_HPP
#define _COMPBIND_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file compBind.hpp

#include <metaProgramming/refOrVal.hpp>
#include <tensor/component.hpp>
#include <tensor/expr.hpp>

namespace nissa
{
  template <typename E,
	    typename BC,
	    bool EvalToConstFlag>
  struct _CompBinderFactory;
  
  template <typename E,
	    typename...TC,
	    typename F,
	    ExprFlags _Flags,
	    typename... BCs,
	    bool EvalToConstFlag>
  struct _CompBinderFactory<Expr<E,TensorComps<TC...>,F,_Flags> // E
			    ,TensorComps<BCs...>,              // BC
			    EvalToConstFlag>
  {
    using BoundExpressionComponents=
      TensorComps<TC...>;
    
    using BoundExpression=
      ConditionalRef<getStoreByRef<_Flags>,E>;
    
    using _BoundComponents=
      TensorComps<BCs...>;
    
    //using CompBinder CompBinder<BOUND_EXPRESSION,BOUND_COMPONENTS>
    using FilteredComponents=
      TupleFilterOut<_BoundComponents,
		     BoundExpressionComponents>;
    
    static constexpr ExprFlags Flags=
      unsetStoreByRef<addEvalToConstIf<EvalToConstFlag,_Flags>>;
    
    struct _CompBinder :
      Expr<_CompBinder,FilteredComponents,F,Flags>
    {
      using Comps=FilteredComponents;
      
      using BoundComponents=
	_BoundComponents;
      
      using Fund=F;
      
      BoundExpression boundExpression;
      
      BoundComponents boundComponents;
      
      template <typename...TD>
      CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
      decltype(auto) eval(const TD&...td) const
      {
	return boundExpression.eval(std::get<BCs>(boundComponents)...,td...);
      }
      
      PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(eval,CUDA_HOST_DEVICE);
      
      _CompBinder(BoundExpression boundExpression,
		 const BoundComponents& boundComponents) :
	boundExpression(boundExpression),
	boundComponents(boundComponents)
	{
	}
    };
  };
  
  template <typename E,
	    typename BC,
	    bool EvalToConstFlag>
  using CompBinder=
    typename _CompBinderFactory<E,BC,EvalToConstFlag>::_CompBinder;
  
  template <typename E,
	    typename...BCs>
  auto compBind(ExprFeat<E>& eFeat,
		const TensorComps<BCs...>& bc)
  {
    constexpr bool doNotEvalToConst=false;
    
    return CompBinder<E,TensorComps<BCs...>,doNotEvalToConst>(eFeat.deFeat().crtp(),bc);
  }
  
//     :
//   fare struttura più grossa dove definire ordinatamente tutto questo,
//       e type dentro
//           ?
// #define BOUND_EXPRESSION_COMPONENTS 


//           template <typename E, typename... TC, typename F, typename... BCs>
//           struct COMP_BINDER
//       : Expr<CompBinder<E, TensorComps<BCs...>>, typename E::Comps>,
//       typename E::Fund > {
//     using Comps=TupleFilterOut<TensorComps<BCs...>,typename E::Comps>;
    
//     using Fund=typename E::Fund;
//   };
}

#endif
