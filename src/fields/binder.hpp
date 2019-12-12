#ifndef _BINDER_HPP
#define _BINDER_HPP

#include "features.hpp"
#include "tensor.hpp"

namespace nissa
{
  /// Binder a component or more than one
  template <typename T,    // Type of the reference to bind
	    typename...C>  // Type of the components to bind
  struct Binder : public Subscribable<Binder<T,C...>>
  {
    /// Reference to bind
    ref_or_val_t<T> ref;
    
    /// Components to bind
    TensComps<ref_or_val_t<C>...> vals;
    
    /// Access to the reference passing all bound components, and more
    template <typename...Tail>
    decltype(auto) eval(Tail&&...tail) const
    {
      return ref(std::get<C>(vals)...,std::forward<Tail>(tail)...);
    }
    
    /// Call a function using the reference, the bound components, and all passed ones
    template <typename F,       // Type of the function
	      typename...Tail>  // Type of the other components
    decltype(auto) call(F&& f,                    ///< Function to call
			Tail&&...tail) const      ///< Other components to pass
    {
      return f(ref,std::get<C>(vals)...,std::forward<Tail>(tail)...);
    }
    
    /// Construct the binder from a reference and components
    Binder(T&& ref,     ///< Reference
	   C&&...vals)  ///< Components
      : ref(ref),vals{vals...}
    {
    }
  };
  
  /// Creates a binder, using the passed reference and component
  template <typename T,          // Reference type
	    typename...Tp>       // Components type
  auto bindComp(T&& ref,       ///< Reference
		Tp&&...comps)  ///< Components
  {
    return Binder<T,Tp...>(std::forward<T>(ref),std::forward<Tp>(comps)...);
  }
  
  // /// Creates a binder, using the passed reference and component
  // template <typename T,    // Reference type
  // 	    typename...Tp> // Components type
  // auto bind2(T&& ref,       ///< Reference
  // 	     Tp&&...comps)  ///< Components
  // {
  //   return Binder<T,Tp...>(std::forward<T>(ref),std::forward<Tp>(comps)...);
  // }
  
  /// Properties of the components of a binder
  template <typename T,
	    typename...C>
  struct CompsTraits<Binder<T,C...>>
  {
    /// Actual type
    using Type=Binder<T,C...>;
    
    /// Type of the bound object
    using Ref=T;
    
    /// Bound components
    using BoundComps=std::tuple<std::remove_reference_t<C>...>;
    
    /// Components of the reference
    using RefComps=typename CompsTraits<typename std::remove_cv_t<std::remove_reference_t<T>>>::Comps;
    
    /// Components
    using Comps=TupleFilterOut<BoundComps,RefComps>;
  };
}

#endif
