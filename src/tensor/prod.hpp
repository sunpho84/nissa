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
    
    return IndependentComponents<Comps1,Comps2>();//Transp<typename CH::E,CH::Flags>(std::forward<_E>(e));
  }
  
}

#endif
