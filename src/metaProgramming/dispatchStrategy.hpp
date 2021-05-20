#ifndef _DISPATCH_STRATEGY_HPP
#define _DISPATCH_STRATEGY_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file disaptchStrategy.hpp

namespace nissa
{
  /// Declares a strategy to be used for tag dispatch, see example on loopOnAllComponents
#define DECLARE_DISPATCH_STRATEGY(NAME,FLAG_NAME)	\
  using NAME=						\
    std::integral_constant<int,FLAG_NAME>*
}

#endif
