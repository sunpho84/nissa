#ifndef _DAG_HPP
#define _DAG_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file dag.hpp

#include <metaProgramming/universalReference.hpp>
#include <tensor/expr.hpp>
#include <tensor/refCatcher.hpp>

namespace nissa
{
  template <typename _E>
  auto dag(_E&& e,
	   UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR)
  {
    return transp(conj(std::forward<_E>(e)));
  }
}

#endif
