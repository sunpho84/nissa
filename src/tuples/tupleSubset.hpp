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
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
  auto tupleGetMany(const std::tuple<Tin...>& in)
  {
    return std::make_tuple(std::get<Tout>(in)...);
  }
  
  /// Fill a tuple with the elements from another tuple
  ///
  /// Note that all elements of out must be in in, but the opposite is
  /// not necessary
  template <typename...Tout,
	    typename...Tin>
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
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
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
  auto tupleGetSubset(const std::tuple<Tin...>& in)
  {
    TpOut out;
    
    tupleFillWithSubset(out,in);
    
    return out;
  }
  
  /////////////////////////////////////////////////////////////////
  
  /// Fill a tuple with the elements from another tuple
  ///
  /// Note that all elements of in must be in out, but the opposite is
  /// not necessary
  template <typename...Tout,
	    typename...Tin>
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
  void tupleSubsetFillWith(std::tuple<Tout...>& out,
			   const std::tuple<Tin...>& in)
  {
    ((std::get<Tin>(out)=std::get<Tin>(in)),...);
  }
}

#endif
