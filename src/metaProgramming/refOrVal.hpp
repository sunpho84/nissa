#ifndef _REFORVAL_HPP
#define _REFORVAL_HPP

#include <type_traits>

namespace nissa
{
  /// If the type is an l-value reference, provide the type T&, otherwise wih T
  template <typename T>
  using ref_or_val_t=std::conditional_t<std::is_lvalue_reference<T>::value,T&,T>;
}

#endif
