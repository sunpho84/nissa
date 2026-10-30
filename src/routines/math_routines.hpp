#ifndef _MATH_ROUTINES_HPP
#define _MATH_ROUTINES_HPP

#include <algorithm>
#include <functional>
#include <concepts>
#include <limits>
#include <cstdint>
#include <array>
#include <utility>

#include <base/debug.hpp>

#include "new_types/complex.hpp"

//Pi
#ifndef M_PI
# define M_PI           3.14159265358979323846
#endif
//sqrt(2)
#define RAD2 1.414213562373095048801688724209l

namespace nissa
{
  /// Low [0] and high [1] 32-bit halves, each held in a uint64_t for arithmetic
  struct Uint64Halves :
    std::array<uint64_t,2>
  {
    /// Split a 64-bit value independently of memory byte order
    CUDA_HOST_AND_DEVICE
    constexpr explicit Uint64Halves(const uint64_t value=0) :
      std::array<uint64_t,2>{value&UINT32_MAX,value>>32}
    {
    }
    
    /// Construct from low and high halves, each of which must fit in 32 bits
    CUDA_HOST_AND_DEVICE
    constexpr Uint64Halves(const uint64_t low,
                          const uint64_t high) :
      std::array<uint64_t,2>{low,high}
    {
    }
    
    /// Reassemble the two valid halves into a 64-bit value
    CUDA_HOST_AND_DEVICE
    constexpr explicit operator uint64_t() const
    {
      return ((*this)[1]<<32)|(*this)[0];
    }
  };
  
  /// Return {low sum,carry}, where carry is zero or one
  CUDA_HOST_AND_DEVICE
  constexpr std::pair<uint64_t,uint64_t> addWithCarry(const uint64_t a,
                                                     const uint64_t b)
  {
    const uint64_t low=a+b;
    
    return {low,low<a};
  }
  
  /// Return {low product,high product}, retaining all 128 product bits
  CUDA_HOST_AND_DEVICE
  constexpr std::pair<uint64_t,uint64_t> multiplyWithCarry(const uint64_t a,
                                                          const uint64_t b)
  {
    const Uint64Halves aHalves{a};
    const Uint64Halves bHalves{b};
    const uint64_t lowProduct=aHalves[0]*bHalves[0];
    
    // The mixed products form a 65-bit column whose carry must be retained
    const auto [cross,crossCarry]=
      addWithCarry(aHalves[0]*bHalves[1],aHalves[1]*bHalves[0]);
    const auto [middle,middleCarry]=addWithCarry(cross,lowProduct>>32);
    const uint64_t high=
      aHalves[1]*bHalves[1]+(middle>>32)+((crossCarry+middleCarry)<<32);
    const Uint64Halves lowHalves{lowProduct&UINT32_MAX,middle&UINT32_MAX};
    
    return {static_cast<uint64_t>(lowHalves),high};
  }
  
  /// Add unsigned integers, rejecting overflow before evaluating the sum
  template <std::unsigned_integral T>
  CUDA_HOST_AND_DEVICE
  constexpr T checkedAdd(const T a,
                         const T b,
                         const char* error="Unsigned addition overflow")
  {
    if(b>std::numeric_limits<T>::max()-a)
      CRASH("%s",error);
    
    return a+b;
  }
  
  /// Multiply unsigned integers, rejecting overflow before evaluating the product
  template <std::unsigned_integral T>
  CUDA_HOST_AND_DEVICE
  constexpr T checkedMultiply(const T a,
                              const T b,
                              const char* error="Unsigned multiplication overflow")
  {
    if(a!=0 and b>std::numeric_limits<T>::max()/a)
      CRASH("%s",error);
    
    return a*b;
  }
  
  double lfact(double n);
  double metro_tresh(double arg);
  int metro_test(double arg);
  
  std::vector<int> factorize(int N);
  
  int log2N(int N);
  CUDA_HOST_AND_DEVICE void matrix_determinant(complex d,complex *m,int n);
  int bitrev(int in,int l2n);
  int find_max_pow2(int a);
  
  //return a bit
  inline bool get_bit(int i,int ibit)
  {return (i>>ibit)&1;}
  
  template <class T>
  T summ(const T& a,
	 const T& b)
  {
    return a+b;
  }
  
  template <class T>
  T nissa_max(const T& a,
	      const T& b)
  {
    return std::max(a,b);
  }
  
  template <class T1,class T2>
  CUDA_HOST_AND_DEVICE
  auto nissa_min(const T1& a,
		 const T2& b)
  {
    return (a<b)?a:b;
  }
  
  template <typename T>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  T sqr(const T& a)
  {
    return a*a;
  }
  
  template <class T>
  void ave_dev(T &ave,T &dev,const T *v,const int n)
  {
    ave=dev=0;
    for(int i=0;i<n;i++)
      {
	ave+=v[i];
	dev+=sqr(v[i]);
      }
    
    ave/=n;
    dev/=n;
    dev-=ave*ave;
    dev=sqrt(dev*n/(n-1));
  }
  
  /// Combine the the passed list of values
  template <typename F,
	    typename T,
	    typename...Ts>
  constexpr T binaryCombine(F&& f,
			    const T& init,
			    Ts&&...list)
  {
    /// Result
    T out=init;
    
    const T l[]{list...};
    
    for(auto i : l)
      out=f(out,i);
    
    return out;
  }
  
  ///Product of the arguments
  template <typename T,
	    typename...Ts>
  constexpr auto productAll(Ts&&...t)
  {
    return binaryCombine(std::multiplies<>(),T{1},std::forward<Ts>(t)...);
  }
  
  ///Sum of the arguments
  template <typename T,
	    typename...Ts>
  constexpr auto sumAll(Ts&&...t)
  {
    return binaryCombine(std::plus<>(),T{0},std::forward<Ts>(t)...);
  }
}

#endif
