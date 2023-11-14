#ifndef _EXPONENTIATOR_HPP
#define _EXPONENTIATOR_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/exponentiator.hpp

#include <expr/cWiseCombine.hpp>

namespace nissa
{
  /// Exponentiator in the case of a complex exponent
  struct ComplExponentiator
  {
    template <typename RealPart,
	      typename ImagPart,
	      typename ReImDiscr>
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    static auto compute(RealPart&& r,
			ImagPart&& i,
			ReImDiscr&& rid)
    {
      if(rid==0)
	return exp(r)*cos(i);
      else
	return exp(r)*sin(i);
    }
  };
  
  /// Exponentiator in the case of a non-complex exponent
  struct Exponentiator
  {
    template <typename E>
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    static auto compute(E&& e)
    {
      return exp(std::forward<E>(e));
    }
  };
  
  /// Catch exp(node)
  template <DerivedFromNode E>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  auto exp(E&& e)
  {
    /// Detect complex
    constexpr bool hasComplComp=
      tupleHasType<typename std::decay_t<E>::Comps,ComplId>;
    
    if constexpr(hasComplComp)
      return
	cWiseCombine<ComplExponentiator>(bindComps(e,std::make_tuple(reIm(0))),
					 bindComps(e,std::make_tuple(reIm(1))),
					 I);
    else
      return
	cWiseCombine<Exponentiator>(std::forward<E>(e));
  }
}

#endif
