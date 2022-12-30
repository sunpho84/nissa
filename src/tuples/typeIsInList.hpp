#ifndef _TYPEISINLIST_HPP
#define _TYPEISINLIST_HPP

#include <tuple>

namespace nissa
{
  /// Predicate returning whether the type is present in the list
  ///
  /// Forward definition
  template <int N,
	    typename Tp>
  struct TypeIsInList;
  
  /// Predicate returning whether the type is present in the list
  template <int N,
	    typename...Tp>
  struct TypeIsInList<N,std::tuple<Tp...>>
  {
    /// Internal implementation
    template <typename T>
    struct t
    {
      /// Predicate result
      static constexpr bool value=((std::is_same<T,Tp>::value+...)==N);
    };
  };
}

#endif
