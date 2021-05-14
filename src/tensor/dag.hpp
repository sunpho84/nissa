#ifndef _DAG_HPP
#define _DAG_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file dag.hpp

#include <utility>

#include <metaProgramming/universalReference.hpp>
#include <tensor/conj.hpp>
#include <tensor/transp.hpp>

namespace nissa
{
  /// Take the hermitean conjugate of e
  template <typename _E,
	    UNPRIORITIZE_DEFAULT_VERSION_TEMPLATE_PARS>
  auto dag(_E&& e,
	   UNPRIORITIZE_DEFAULT_VERSION_ARGS)
  {
    UNPRIORITIZE_DEFAULT_VERSION_ARGS_CHECK;
    
    return
      transp(conj(std::forward<_E>(e)));
  }
  
  /// Take the hermitean conjugate of e when e is a conjugator
  template <typename _E,
	    ENABLE_THIS_TEMPLATE_IF(isConjugator<_E>)>
  auto dag(_E&& e)
  {
    return
      transp(e.nestedExpr);
  }
  
  /// Take the hermitean conjugate of e when e is a transposer
  template <typename _E,
	    ENABLE_THIS_TEMPLATE_IF(isTransposer<_E>)>
  auto dag(_E&& e)
  {
    return
      conj(e.nestedExpr);
  }
}

#endif
