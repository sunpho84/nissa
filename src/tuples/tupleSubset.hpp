#ifndef _TUPLESUBSET_HPP
#define _TUPLESUBSET_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file tupleSubset.hpp

#include <tuple>

#include <metaprogramming/inline.hpp>

namespace nissa
{
  /// Get the list elements from a tuple
  template <typename...Tout,
	    typename...Tin>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  auto tupleGetMany(const std::tuple<Tin...>& in)
  {
    return std::make_tuple(std::get<Tout>(in)...);
  }
  
  /// Fill a tuple with the elements from another tuple
  template <typename...Tout,
	    typename...Tin>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  void tupleFillWithSubset(std::tuple<Tout...>& out,
			   const std::tuple<Tin...>& in)
  {
    ((std::get<Tout>(out)=std::get<Tout>(in)),...);
  }
  
  /// Get the list elements of the passed tuple
  ///
  /// \example
  /// auto tupleGetSubset<std::tuple<int>>(std::make_tuple(1,10.0));
  template <typename TpOut,
	    typename...Tin>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  auto tupleGetSubset(const std::tuple<Tin...>& in)
  {
    TpOut out;
    
    tupleFillWithSubset(out,in);
    
    return out;
  }
}

#endif
