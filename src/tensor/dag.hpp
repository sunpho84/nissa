#ifndef _DAG_HPP
#define _DAG_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file dag.hpp

#include <utility>

#include <metaProgramming/universalReference.hpp>

namespace nissa
{
  /// Take the hermitean conjugate of e
  template <typename _E>
  auto dag(_E&& e,
	   UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR)
  {
    return
      transp(conj(std::forward<_E>(e)));
  }
}

#endif
