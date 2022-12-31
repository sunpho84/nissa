#ifndef _PRODUCER_DECLARATION_HPP
#define _PRODUCER_DECLARATION_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/nodes/producerDeclaration.hpp

#include <utility>

#include <metaprogramming/detectableAs.hpp>

namespace nissa
{
  PROVIDE_DETECTABLE_AS(Producer);
  
  /// Producer
  ///
  /// Forward declaration to capture the components
  template <typename _Ccs,
	    typename _E,
	    typename _Comps,
	    typename _Fund,
	    typename _Is=std::integer_sequence<int,0,1>>
  struct Producer;
}

#endif
