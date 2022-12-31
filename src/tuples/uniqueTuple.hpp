#ifndef _UNIQUETUPLE_HPP
#define _UNIQUETUPLE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file uniqueTuple.hpp

#include <tuple>

namespace nissa
{
  namespace impl
  {
    /// Tuple with unique types from a list
    ///
    /// Forward definition
    template <typename U,
	      typename...Os>
    struct _UniqueTuple;
    
    /// Tuple with unique types from a list
    ///
    /// End case with no type to add
    template <typename...Us>
    struct _UniqueTuple<std::tuple<Us...>>
    {
      /// Provide the list passed as input
      using type=
	std::tuple<Us...>;
    };
    
    /// Tuple with unique types from a list
    ///
    /// Default case processing type O
    template <typename...Us,
	      typename O,
	      typename...Os>
    struct _UniqueTuple<std::tuple<Us...>,
			O,Os...>
    {
      /// Check if the type O is unique
      static constexpr bool isOUnique=
	((not std::is_same_v<Us,O>)&...&true);
      
      /// Compute next list of unique type which might include O
      using NU=
	std::conditional_t<isOUnique,
			   std::tuple<Us...,O>,
			   std::tuple<Us...>>;
      
      /// Get the full list from recursive analysis
      using type=
	typename _UniqueTuple<NU,Os...>::type;
    };
  }
  
  /// Tuple with unique types from a list
  ///
  /// Based on https://stackoverflow.com/a/57528226
  template <typename...Ts>
  using UniqueTuple=
    typename impl::_UniqueTuple<std::tuple<>,Ts...>::type;
  
}

#endif
