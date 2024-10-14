#ifndef _MATH_ROUTINES_HPP
#define _MATH_ROUTINES_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <vector>

#include <metaprogramming/inline.hpp>

namespace nissa
{
  /// Difference with next multiple of N
  template <auto N,
	    typename T>
  constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
  T diffWithNextMultipleOf(const T& x)
  {
    return (N-x%N)%N;
  }
  
  /// Ceil to next multiple of eight
  template <auto N,
	    typename T>
  constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
  T ceilToNextMultipleOf(const T& x)
  {
    return x+diffWithNextMultipleOf<N>(x);
  }
  
  /// Factorize a number
  template <typename T>
  INLINE_FUNCTION
  std::vector<T> factorize(T in)
  {
    std::vector<T> out;
    
    for(T fatt=2;in>1;)
      {
	const T div=in/fatt;
	const T res=in-div*fatt;
	
	if(res!=0)
	  fatt++;
	else
	  {
	    in=div;
	    out.push_back(fatt);
	  }
      }
    
    if(out.size()==0)
      out.push_back(1);
    
    return out;
  }
  
  /// Returns max! / min!
  template <typename T>
  constexpr auto partialFact(const T& min,
			     const T& max)
  {
    /// Result to be returned
    
    T res=1;
    for(T i=min;i<=max;i++)
      res*=i;
    
    return res;
  }
  
  /// Returns N!
  template <typename T>
  constexpr auto fact(const T& N)
  {
    return partialFact((T)1,N);
  }
}

#endif
