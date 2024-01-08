#ifndef _DAGGER_HPP
#define _DAGGER_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/dagger.hpp

#include <expr/conj.hpp>
#include <expr/node.hpp>
#include <expr/transp.hpp>

namespace nissa
{
  /// Take the dagger of an expression
  template <DerivedFromNode _E>
  HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
  decltype(auto) dag(_E&& e)
  {
    using E=std::decay_t<_E>;
    
    if constexpr(isConjugator<E>)
      return transpose(FORWARD_MEMBER_VAR(_E,std::forward<E>(e),template subNode<0>));
    else
      if constexpr(isTransposer<E>)
	return conj(FORWARD_MEMBER_VAR(_E,std::forward<E>(e),transpExpr));
      else
	return transp(conj(std::forward<_E>(e)));
  }
}

#endif
