#ifndef _PROD_HPP
#define _PROD_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file prod.hpp

#include <metaProgramming/universalReference.hpp>
#include <tensor/unaryExpr.hpp>

namespace nissa
{
  namespace internal
  {
    template <typename TT,
	      typename TP1,
	      typename TP2>
    struct _ProdComponents;
    
    template <typename...S,
	      RwCl...RCA,
	      int...Which,
	      typename TP1,
	      typename TP2>
    struct _ProdComponents<TensorComps<TensorComp<S,RCA,Which>...>,
			   TP1,
			   TP2>
    {
      template <typename  _S,
		int _Which>
      static constexpr bool contract=
	tupleHasType<TP1,TensorComp<_S,CLN,_Which>,1> and
	tupleHasType<TP2,TensorComp<_S,ROW,_Which>,1>;
      
      template <typename TC>
      struct _Res;
      
      template <typename _S,
		int _Which>
      struct _Res<TensorComp<_S,ROW,_Which>>
      {
	using ThisRow=
	  TensorComp<_S,ROW,_Which>;
	
	using ThisCln=
	  TensorComp<_S,CLN,_Which>;
	
	static constexpr bool resHasRow=
	  tupleHasType<TP1,ThisRow,1> or (tupleHasType<TP2,ThisRow,1> and not contract<_S,_Which>);
	
	static constexpr bool resHasCln=
	  tupleHasType<TP2,ThisCln,1> or (tupleHasType<TP1,ThisCln,1> and not contract<_S,_Which>);
	
	using type=
	  TupleCat<std::conditional_t<resHasRow,TensorComps<ThisRow>,TensorComps<>>,
		   std::conditional_t<resHasCln,TensorComps<ThisCln>,TensorComps<>>>;
      };
      
      template <typename _S,
		int _Which>
      struct _Res<TensorComp<_S,ANY,_Which>>
      {
	using type=
	  TensorComps<TensorComp<_S,ANY,_Which>>;
      };
      
      using type=
	TupleCat<typename _Res<TensorComp<S,RCA,Which>>::type...>;
    };
  }
  
  template <typename TC1,
	    typename TC2>
  using ProdComps=
    typename internal::_ProdComponents<IndependentComponents<TC1,TC2>,
				       TC1,TC2>::type;
  
  /////////////////////////////////////////////////////////////////
  
  /// Product of two expressions
  template <typename T1,
	    typename T2,
	    typename _Comps,
	    typename _Fund>
  struct Prod : Expr<Prod<T1,T2,_Comps,_Fund>,
		     _Comps,
		     _Fund>
  {
    using Comps=
      _Comps;
        
    using Fund=
      _Fund;
  };
  
  template <typename _E1,
	    typename _E2>
  auto prod(_E1&& e1,
	    _E2&& e2,
	    UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR)
  {
    // using CH1=
    //   RefCatcherHelper<_E1,decltype(e1)>;
    
    // using CH2=
    //   RefCatcherHelper<_E2,decltype(e2)>;
    
    // using F1=
    //   typename CH1::E::Fund;
    
    // using F2=
    //   typename CH2::E::Fund;
    
    // using F=
    //   decltype(F1()*F2());
    
    // using E1=
    //   typename CH1::E;
    
    // using E2=
    //   typename CH2::E;
    
    // return ProdComps<typename E1::Comps,typename E2::Comps>();//Prod<E1,E2,ProdComps<E1,E2>, ExprFlags _Flags>
    /////////////////////////////////////////////////////////////////
  }
}

#endif
