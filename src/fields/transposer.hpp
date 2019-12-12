#ifndef _TRANSPOSER_HPP
#define _TRANSPOSER_HPP

#include <fields/features.hpp>
#include <base/metaprogramming.hpp>

namespace nissa
{
  /// Transpose the reference
  template <typename T>    // Type of the reference to transpose
  struct Transposer : public Subscribable<Transposer<T>>
  {
    /// Reference to bind
    ref_or_val_t<T> ref;
    
    /// Access to the reference transposeing all passed values
    template <typename...Cp>
    decltype(auto) eval(Cp&&...c) const
    {
      return ref(c.transp()...);
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD(eval);
    
    /// Construct the transposer a reference
    Transposer(T&& ref)     ///< Reference
      : ref(ref)
    {
    }
  };
  
  /// Creates a transposer, using the passed reference
  template <typename T>                  // Reference type
  Transposer<T> transpose(T&& ref)      ///< Reference
  {
    return std::forward<T>(ref);
  }
  
  /// Transposed components
  ///
  ///  Forward implementation
  template <typename T>
  struct _TransposedComps;
  
  /// Transposed components
  template <typename...T>
  struct _TransposedComps<std::tuple<T...>>
  {
    /// Input components
    using Comps=std::tuple<T...>;
    
    /// Result
    using type=std::tuple<std::conditional_t<TypeIsInList<1,Comps>::template t<typename T::Transp>::value,T,typename T::Transp>...>;
  };
  
  /// Properties of components of transposer
  template <typename T>
  struct CompsTraits<Transposer<T>>
  {
    /// Actual type
    using Type=Transposer<T>;
    
    /// Reference type
    using Ref=T;
    
    /// Components of the reference
    using RefComps=typename CompsTraits<std::remove_reference_t<T>>::Comps;
    
    /// Transposer components
    using Comps=typename _TransposedComps<RefComps>::type;
  };
}

#endif
