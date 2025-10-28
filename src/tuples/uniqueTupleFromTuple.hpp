#ifndef _UNIQUETUPLEFROMTUPLE_HPP
#define _UNIQUETUPLEFROMTUPLE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file uniqueTupleFromTuple.hpp

#include <tuples/uniqueTuple.hpp>

namespace nissa
{
  namespace impl
  {
    /// Returns a tuple containing only the unique elements from a tuple
    ///
    /// Internal implementation, forward definition
    template <typename Tp>
    struct _UniqueTupleFromTuple;
    
    /// Returns a tuple containing only the unique elements from a tuple
    ///
    /// Internal implementation
    template <typename...T>
    struct _UniqueTupleFromTuple<std::tuple<T...>>
    {
      using type=UniqueTuple<T...>;
    };
  }
  
  /// Returns a tuple containing only the unique elements from a tuple
  ///
  /// Based on https://stackoverflow.com/a/57528226
  template <typename Tp>
  using UniqueTupleFromTuple=
    typename impl::_UniqueTupleFromTuple<Tp>::type;
}

#endif
