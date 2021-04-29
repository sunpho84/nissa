#ifndef _REFORVAL_HPP
#define _REFORVAL_HPP

#include <type_traits>

namespace nissa
{
  namespace details
  {
    template <bool IsRef,
	      typename T>
    using _ConditionalRef=
      std::conditional_t<IsRef,T&,T>;
  }
  
  template <bool IsRef,
	    typename T>
  using ConditionalRef=
    details::_ConditionalRef<IsRef,std::remove_reference_t<T>>;
}

#endif
