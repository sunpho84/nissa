#ifndef _METAPROGRAMMING_HPP
#define _METAPROGRAMMING_HPP

#include <type_traits>
#include <utility>

namespace nissa
{
  
#if __cplusplus > 201402L
  
  using std::as_const;
  
#else
  
  template <typename T>
  constexpr const T& as_const(T& t) noexcept
  {
    return t;
  }
  
#endif
  
  /// Returns true if T is a const lvalue reference
  template <typename T>
  constexpr bool is_const_lvalue_reference_v=std::is_lvalue_reference<T>::value and std::is_const<std::remove_reference_t<T>>::value;
  
  /// Returns the type without "const" attribute if it is a reference
  template <typename T>
  decltype(auto) remove_const_if_ref(T&& t)
  {
    using Tv=std::remove_const_t<std::remove_reference_t<T>>;
    
    return (std::conditional_t<is_const_lvalue_reference_v<T>,Tv&,Tv>)t;
  }
  
  /// If the type is an l-value reference, provide the type T&, otherwise wih T
  template <typename T>
  using ref_or_val_t=std::conditional_t<std::is_lvalue_reference<T>::value,T&,T>;
  
  /// Provides also a non-const version of the method \c NAME
  ///
  /// See
  /// https://stackoverflow.com/questions/123758/how-do-i-remove-code-duplication-between-similar-const-and-non-const-member-func
  /// A const method NAME must be already present Example
  ///
  /// \code
  // class ciccio
  /// {
  ///   double e{0};
  ///
  /// public:
  ///
  ///   const double& get() const
  ///   {
  ///     return e;
  ///   }
  ///
  ///   PROVIDE_ALSO_NON_CONST_METHOD(get);
  /// };
  /// \endcode
#define PROVIDE_ALSO_NON_CONST_METHOD(NAME)				\
  /*! Overload the \c NAME const method passing all args             */ \
  template <typename...Ts> /* Type of all arguments                  */	\
  decltype(auto) NAME(Ts&&...ts) /*!< Arguments                      */ \
  {									\
    return remove_const_if_ref(as_const(*this).NAME(std::forward<Ts>(ts)...)); \
  }
  
  /// Implements the CRTP pattern
  template <typename T>
  struct Crtp
  {
    /// Crtp access the type
    const T& crtp() const
    {
      return *static_cast<const T*>(this);
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD(crtp);
  };
}

#endif
