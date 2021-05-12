#ifndef _PROD_HPP
#define _PROD_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file prod.hpp

#include <metaProgramming/universalReference.hpp>
#include <tensor/expr.hpp>
#include <tensor/refCatcher.hpp>

namespace nissa
{
  /// Product of two expressions
  template <typename T1,
	    typename T2,
	    ExprFlags _Flags>
  struct Prod : Expr<Prod<T1,T2,_Flags>,
		     typename T1::Comps,
		     typename T1::Fund,
		     unsetEvalToRef<setStoreByRefTo<false,_Flags>>>
  {
  };
  
  namespace internal
  {
    template <typename TH,
	      typename TT,
	      typename TP1,
	      typename TP2>
    struct _ProdComponents;
    
    template <typename TH,
	      typename TP1,
	      typename TP2>
    struct _ProdComponents<TH,TensorComps<>,TP1,TP2>
    {
      using type=TH;
    };
    
    template <typename...TcHead,
	      typename ThisS,
	      int ThisWhich,
	      typename...TcTail,
	      typename TP1,
	      typename TP2>
    struct _ProdComponents<TensorComps<TcHead...>,
			   TensorComps<TensorComp<ThisS,ANY,ThisWhich>,
				       TcTail...>,
			   TP1,
			   TP2>
    {
      using type=
	typename _ProdComponents<TensorComps<TcHead...,TensorComp<ThisS,ANY,ThisWhich>>,
				 TensorComps<TcTail...>,
				 TP1,
				 TP2>::type;
    };
    
    template <typename Tc,
	      typename TP>
    struct TensorCompsContainsRowClnComp;
    
    template <typename S,
	      int Which,
	      typename...TP>
    struct TensorCompsContainsRowClnComp<TensorComp<S,ROW,Which>,TensorComps<TP...>>
    {
      template <RwCl RC>
      static constexpr bool hasOne=
	((std::is_same_v<TensorComp<S,RC,Which>,TP>)+...)==1;
      
      static constexpr bool hasRow=
	hasOne<ROW>;
      
      static constexpr bool hasCln=
	hasOne<CLN>;
      
      static constexpr bool hasRowAndCln=
	hasRow and hasCln;
    };
    
    template <typename Tc,
	      typename TP>
    constexpr bool hasRowComp=
      TensorCompsContainsRowClnComp<Tc,TP>::hasRow;
    
    template <typename Tc,
	      typename TP>
    constexpr bool hasClnComp=
      TensorCompsContainsRowClnComp<Tc,TP>::hasCln;
    
    template <typename Tc,
	      typename TP>
    constexpr bool hasRowAndClnComp=
      TensorCompsContainsRowClnComp<Tc,TP>::hasRowAndCol;
    
    template <typename...TcHead,
	      typename ThisS,
	      int ThisWhich,
	      typename...TcTail,
	      typename TP1,
	      typename TP2>
    struct _ProdComponents<TensorComps<TcHead...>,
			   TensorComps<TensorComp<ThisS,ROW,ThisWhich>,
				       TcTail...>,
			   TP1,
			   TP2>
    {
      using ThisComp=
	TensorComp<ThisS,ROW,ThisWhich>;
      
      static constexpr bool contract=
	hasClnComp<ThisComp,TP1> and
	hasRowComp<ThisComp,TP2>;
      
      template <RwCl RC>
      struct _ResHas;
      
      template <>
      struct _ResHas<ROW>
      {
	static constexpr bool value=
	  hasRowComp<ThisComp,TP1> or (hasRowComp<ThisComp,TP2> and not contract);
      };
      
      template <>
      struct _ResHas<CLN>
      {
	static constexpr bool value=
	  hasClnComp<ThisComp,TP2> or (hasClnComp<ThisComp,TP1> and not contract);
      };
      
      template <RwCl RC>
      static constexpr bool resHas=
	_ResHas<RC>::value;
      
      template <RwCl RC>
      using _ThisCompConditionallyInclude=
	std::conditional_t<resHas<RC>,TensorComps<TensorComp<ThisS,RC,ThisWhich>>,TensorComps<>>;
      
      using type=
	typename _ProdComponents<TupleCat<TensorComps<TcHead...>,
					  _ThisCompConditionallyInclude<ROW>,
					  _ThisCompConditionallyInclude<CLN>>,
				 TensorComps<TcTail...>,
				 TP1,
				 TP2>::type;
    };
  }
  
  template <typename TC1,
	    typename TC2>
  using ProdComps=
    typename internal::_ProdComponents<TensorComps<>,
				       IndependentComponents<TC1,TC2>,
				       TC1,TC2>::type;

  template <typename _E1,
	    typename _E2>
  auto prod(_E1&& e1,
	    _E2&& e2,
	    UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR)
  {
    using CH1=
      RefCatcherHelper<_E1,decltype(e1)>;
    
    using CH2=
      RefCatcherHelper<_E2,decltype(e2)>;
    
    using F1=
      typename CH1::E::Fund;
    
    using F2=
      typename CH2::E::Fund;
    
    using F=
      decltype(F1()*F2());
    
    using Comps1=
      typename CH1::E::Comps;
    
    using Comps2=
      typename CH2::E::Comps;
    
    return ProdComps<Comps1,Comps2>();
  }
}

#endif
