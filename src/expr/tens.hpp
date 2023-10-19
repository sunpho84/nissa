#ifndef _TENS_HPP
#define _TENS_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/tens.hpp

#include <expr/dynamicTens.hpp>
#include <expr/stackTens.hpp>

namespace nissa
{
  /// Provide a tensor, reverting to stack if possible
  template <typename C,
	    typename Fund>
  using Tens=
    std::conditional_t<compsAreDynamic<C>,DynamicTens<C,Fund>,StackTens<C,Fund>>;
}

#endif
