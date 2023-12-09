#ifndef _MATH_ROUTINES_HPP
#define _MATH_ROUTINES_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <algorithm>
#include <cmath>
#include <functional>

#include <metaprogramming/inline.hpp>

namespace nissa
{
  double lfact(double n);
  double metro_tresh(double arg);
  int factorize(int *list,int N);
  int log2N(int N);
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
  
  template <class T>
  T cube(T a)
  {return a*a*a;};
  
  /// Integral power
  template <class T,
	    std::incrementable I>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  T pow(const T& in,
	const I& n)
  {
    T out=1;
    
    for(I i=0;i<n;i++)
      out*=in;
    
    return out;
  };
  
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
    
    return out;
  }
}

#endif
