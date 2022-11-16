#ifndef _ASSIGN_HPP
#define _ASSIGN_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file assign.hpp
///
/// \brief Implements assignment of expressions

#include <metaProgramming/inliner.hpp>
#include <tensor/expr.hpp>
#include <tensor/loopOnAllComponents.hpp>

namespace nissa
{
  template <typename EL,
	    typename ER>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
  EL& assign(EL& lhs,
	     const ER& rhs)
  {
    loopOnAllComponents<typename ER::Comps>(lhs.getDynamicSizes(),
					    [&](const auto&...comps) INLINE_ATTRIBUTE
					    {
					      lhs(comps...)=rhs(comps...);
					    });
    
    return lhs;
  }
}

#endif
