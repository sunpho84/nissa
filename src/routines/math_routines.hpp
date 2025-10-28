#ifndef _MATH_ROUTINES_HPP
#define _MATH_ROUTINES_HPP

#include <algorithm>
#include <functional>

#include "new_types/complex.hpp"

//Pi
#ifndef M_PI
# define M_PI           3.14159265358979323846
#endif
//sqrt(2)
#define RAD2 1.414213562373095048801688724209l

namespace nissa
{
  double lfact(double n);
  double metro_tresh(double arg);
  int metro_test(double arg);
  int factorize(int *list,int N);
  int log2N(int N);
  CUDA_HOST_AND_DEVICE void matrix_determinant(complex d,complex *m,int n);
  int bitrev(int in,int l2n);
  int find_max_pow2(int a);
  
  //return a bit
  inline bool get_bit(int i,int ibit)
  {return (i>>ibit)&1;}
  
  template <class T>
  T summ(const T& a,const T& b)
  {return a+b;}
  
  template <class T>
  T nissa_max(const T& a,const T& b)
  {
    return std::max(a,b);
  }
  
  template <class T1,class T2>
  CUDA_HOST_AND_DEVICE
  auto nissa_min(const T1& a,const T2& b)
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
