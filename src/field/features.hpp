#ifndef _FEATURES_HPP
#define _FEATURES_HPP

#include <base/metaprogramming.hpp>

namespace nissa
{
  /// Component traits
  ///
  /// Forward definition
  template <typename T>
  struct CompsTraits;
  
  /// Subscribable feature
  template <typename T>
  struct Subscribable : public Crtp<T>
  {
    /// Number of components
    static constexpr int NComps=std::tuple_size<typename CompsTraits<T>::Comps>::value;
    
    /// Single component access via subscribe operator
    template <typename S>                      // Type of the subscribed component
    decltype(auto) operator[](S&& s) const     ///< Subscribed component
    {
      return (*this)(std::forward<S>(s));
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD(operator[]);
    
    /// Access to inner data with any order
    template <typename...Cp,
	      std::enable_if_t<sizeof...(Cp)!=NComps,void*> =nullptr>
    decltype(auto) operator()(Cp&&...comps) const ///< Components
    {
      return bindComp(this->crtp(),std::forward<Cp>(comps)...);
    }
    
    /// Access to inner data with any order
    template <typename...Cp,
	      std::enable_if_t<sizeof...(Cp)==NComps,void*> =nullptr>
    decltype(auto) operator()(Cp&&...comps) const ///< Components
    {
      return this->crtp().eval(std::forward<Cp>(comps)...);
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD(operator());
  };
}

#endif
