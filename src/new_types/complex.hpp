#ifndef _COMPLEX_HPP
#define _COMPLEX_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <stdio.h>

//real/imag
#define RE 0
#define IM 1

#include <base/metaprogramming.hpp>

namespace nissa
{
  typedef double complex[2];
  typedef complex quad_u1[NDIM];
  
  typedef float single_complex[2];
  
  //////////////////////////////////////////////////////////
  
  inline double real_part_of_complex_prod(const complex a,const complex b)
  {return a[0]*b[0]-a[1]*b[1];};
  
  /// Re(a,b)
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  double real_part_of_complex_scalar_prod(const A& a,
					  const B& b)
  {
    return a[0]*b[0]+a[1]*b[1];
  };
  
  //print
  inline void complex_print(const complex a)
  {printf("(%16.16lg,%16.16lg)\n",a[0],a[1]);}
  
  /// a=b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_copy(A&& a,
		    const B& b)
  {
    a[0]=b[0];
    a[1]=b[1];
  }
  
  inline void complex_copy_from_single_complex(complex a,const single_complex b)
  {
    a[0]=b[0];
    a[1]=b[1];
  }
  inline void single_complex_copy_from_complex(single_complex a,const complex b)
  {
    a[0]=b[0];
    a[1]=b[1];
  }
  
  /// a=(b,0)
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_put_to_real(A&& a,
			   const double& b)
  {
    a[0]=b;
    a[1]=0;
  }
  
  inline void complex_put_to_imag(complex a,const double b)
  {
    a[0]=0;
    a[1]=b;
  }
  
  /// a=0
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_put_to_zero(A&& a)
  {
    complex_put_to_real(a,0);
  }
  
  /// a=~b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_conj(A&& a,
		    const B& b)
  {
    a[0]=b[0];
    a[1]=-b[1];
  }
  
  //Assign minus the conj
  inline void complex_minus_conj(complex a,const complex b)
  {
    a[0]=-b[0];
    a[1]=b[1];
  }
  
  /// a=b+c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_summ(A&& a,
		    const B& b,
		    const C& c)
  {
    a[0]=b[0]+c[0];
    a[1]=b[1]+c[1];
  }
  
  /// a=b+i*c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_isumm(A&& a,
		     const B& b,
		     const C& c)
  {
    a[0]=b[0]-c[1];
    a[1]=b[1]+c[0];
  }
  
  /// a+=i*b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_isummassign(A&& a,
			   const B& b)
  {
    complex_isumm(a,a,b);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_summ_conj2(A&& a,
			  const B& b,
			  const C& c)
  {
    a[0]=b[0]+c[0];
    a[1]=b[1]-c[1];
  }
  
  inline void complex_summ_conj1(complex a,const complex b,const complex c)
  {complex_summ_conj2(a,c,b);}
  
  /// a=b-c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_subt(A&& a,
		    const B& b,
		    const C& c)
  {
    a[0]=b[0]-c[0];
    a[1]=b[1]-c[1];
  }
  
  /// a=b-i*c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_isubt(A&& a,
		    const B& b,
		    const C& c)
  {
    a[0]=b[0]+c[1];
    a[1]=b[1]-c[0];
  }
  
  /// a-=i*b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_isubtassign(A&& a,
			   const B& b)
  {
    complex_isubt(a,a,b);
  }
  
  CUDA_HOST_AND_DEVICE inline void complex_subt_conj2(complex a,const complex b,const complex c)
  {
    a[0]=b[0]-c[0];
    a[1]=b[1]+c[1];
  }
  inline void complex_subt_conj1(complex a,const complex b,const complex c)
  {
    a[0]=+b[0]-c[0];
    a[1]=-b[1]-c[1];
  }
  
  /// a+=b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_summassign(A&& a,
			  const B& b)
  {
    complex_summ(a,a,b);
  }
  
  /// a-=b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_subtassign(A&& a,
			  const B& b)
  {
    complex_subt(a,a,b);
  }
  
  //put to exp
  CUDA_HOST_AND_DEVICE inline void complex_iexp(complex out,const double arg)
  {sincos(arg,out+IM,out+RE);}
  
  /// Prod with real
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_prod_double(A&& a,
			   const B& b,
			   const double& c)
  {
    a[RE]=b[RE]*c;
    a[IM]=b[IM]*c;
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_prodassign_double(A&& a,
				 const double& c)
  {
    complex_prod_double(a,a,c);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_prod_idouble(A&& a,
			  const B& b,
			  const C& c)
  {
    const double d=-b[IM]*c;
    
    a[IM]=b[RE]*c;
    a[RE]=d;
  }
  
  CUDA_HOST_AND_DEVICE inline void complex_prodassign_idouble(complex a,const double b) {complex_prod_idouble(a,a,b);}
  
  //summ the prod with real
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_summ_the_prod_double(A&& a,
				    const B& b,
				    const double& c)
  {
    const auto t=b[0]*c;
    a[1]+=b[1]*c;
    a[0]+=t;
  }
  
  CUDA_HOST_AND_DEVICE inline void complex_subt_the_prod_double(complex a,const complex b,const double c)
  {
    const double t=b[0]*c;
    a[1]-=b[1]*c;
    a[0]-=t;
  }
  
  //summ the prod with imag
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_summ_the_prod_idouble(A&& a,
				     const B& b,
				     const double& c)
  {
    const auto t=b[1]*c;
    
    a[1]+=b[0]*c;
    a[0]-=t;
  }
  
  /// a-=b*i*c
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_subt_the_prod_idouble(A&& a,
				     const B& b,
				     const double& c)
  {
    const double t=b[1]*c;
    
    a[1]-=b[0]*c;
    a[0]+=t;
  }
  
  //take a linear combination
  inline void complex_linear_comb(complex a,const complex b,const double cb,const complex c,const double cc)
  {
    a[RE]=b[RE]*cb+c[RE]*cc;
    a[IM]=b[IM]*cb+c[IM]*cc;
  }
  
  /// Summ to the output the product of two complex number
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_summ_the_prod(A&& a,const B& b,const C& c)
  {
    const auto t=b[0]*c[0]-b[1]*c[1];
    
    a[1]+=b[0]*c[1]+b[1]*c[0];
    a[0]+=t;
  }
  
  CUDA_HOST_AND_DEVICE inline void single_complex_summ_the_prod(single_complex a,const single_complex b,const single_complex c)
  {
    const double t=b[0]*c[0]-b[1]*c[1];
    a[1]+=b[0]*c[1]+b[1]*c[0];
    a[0]+=t;
  }
  
  /// Subt from the output the product of two complex number
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_subt_the_prod(A&& a,
			     const B& b,
			     const C& c)
  {
    const auto t=b[0]*c[0]-b[1]*c[1];
    
    a[1]-=b[0]*c[1]+b[1]*c[0];
    a[0]-=t;
  }
  
  /// a+=b*~c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_summ_the_conj2_prod(A&& a,
				   const B& b,
				   const C& c)
  {
    const auto t=+b[0]*c[0]+b[1]*c[1];
    
    a[1]+=-b[0]*c[1]+b[1]*c[0];
    a[0]+=t;
  }
  
  inline void single_complex_summ_the_conj2_prod(single_complex a,const single_complex b,const single_complex c)
  {
    const double t=+b[0]*c[0]+b[1]*c[1];
    a[1]+=-b[0]*c[1]+b[1]*c[0];
    a[0]+=t;
  }
  
  /// a+=~b*c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_summ_the_conj1_prod(A&& a,
				   const B& b,
				   const C& c)
  {
    complex_summ_the_conj2_prod(a,c,b);
  }
  
  inline void single_complex_summ_the_conj1_prod(single_complex a,const single_complex b,const single_complex c)
  {single_complex_summ_the_conj2_prod(a,c,b);}
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_summ_the_conj_conj_prod(A&& a,
				       const B& b,
				       const C& c)
  {
    const auto t=+b[0]*c[0]-b[1]*c[1];
    a[1]+=-b[0]*c[1]-b[1]*c[0];
    a[0]+=t;
  }
  
  /// a-=b*~c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_subt_the_conj2_prod(A&& a,
				   const B& b,
				   const C& c)
  {
    const auto t=+b[0]*c[0]+b[1]*c[1];
    
    a[1]-=-b[0]*c[1]+b[1]*c[0];
    a[0]-=t;
  }
  
  CUDA_HOST_AND_DEVICE inline void single_complex_subt_the_conj2_prod(single_complex a,const single_complex b,const single_complex c)
  {
    const double t=+b[0]*c[0]+b[1]*c[1];
    a[1]-=-b[0]*c[1]+b[1]*c[0];
    a[0]-=t;
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_subt_the_conj1_prod(A&& a,
				   const B& b,
				   const C& c)
  {
    complex_subt_the_conj2_prod(a,c,b);
  }
  
  CUDA_HOST_AND_DEVICE inline void single_complex_subt_the_conj1_prod(single_complex a,const single_complex b,const single_complex c)
  {single_complex_subt_the_conj2_prod(a,c,b);}
  CUDA_HOST_AND_DEVICE inline void complex_subt_the_conj_conj_prod(complex a,const complex b,const complex c)
  {
    const double t=+b[0]*c[0]-b[1]*c[1];
    a[1]-=-b[0]*c[1]-b[1]*c[0];
    a[0]-=t;
  }
  
  /// Product of two complex number
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_complex_prod(A&& a,
			   const B& b,
			   const C& c)
  {
    a[0]=b[0]*c[0]-b[1]*c[1];
    a[1]=b[0]*c[1]+b[1]*c[0];
  }
  
  inline void unsafe_single_single_complex_prod(single_complex a,const single_complex b,const single_complex c)
  {
    a[0]=b[0]*c[0]-b[1]*c[1];
    a[1]=b[0]*c[1]+b[1]*c[0];
  }
  CUDA_HOST_AND_DEVICE inline void unsafe_single_complex_prod(single_complex a,const single_complex b,const single_complex c)
  {
    a[0]=b[0]*c[0]-b[1]*c[1];
    a[1]=b[0]*c[1]+b[1]*c[0];
  }
  
  //Minus the product
  inline void unsafe_complex_prod_minus(complex a,const complex b,const complex c)
  {
    a[0]=-(b[0]*c[0]-b[1]*c[1]);
    a[1]=-(b[0]*c[1]+b[1]*c[0]);
  }
  
  /// a=b*~c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_complex_conj2_prod(A&& a,
				 const B& b,
				 const C& c)
  {
    a[0]=+b[0]*c[0]+b[1]*c[1];
    a[1]=-b[0]*c[1]+b[1]*c[0];
  }
  
  //Minus the former
  inline void unsafe_complex_conj2_prod_minus(complex a,const complex b,const complex c)
  {
    a[0]=-(b[0]*c[0]+b[1]*c[1]);
    a[1]=-(-b[0]*c[1]+b[1]*c[0]);
  }
  
  /// a=~b*c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_complex_conj1_prod(A&& a,const B& b,const C& c)
  {
    unsafe_complex_conj2_prod(a,c,b);
  }
  
  inline void unsafe_complex_conj1_prod_minus(complex a,const complex b,const complex c)
  {unsafe_complex_conj2_prod_minus(a,c,b);}
  
  //The product of the conjugate of two complex numbers
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_complex_conj_conj_prod(A&& a,
				     const B& b,
				     const C& c)
  {
    a[0]=+b[0]*c[0]-b[1]*c[1];
    a[1]=-b[0]*c[1]-b[1]*c[0];
  }
  
  //Minus the former
  inline void unsafe_complex_conj_conj_prod_minus(complex a,const complex b,const complex c)
  {
    a[0]=-b[0]*c[0]+b[1]*c[1];
    a[1]=+b[0]*c[1]+b[1]*c[0];
  }
  
  /// The product of two complex number
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_complex_prod(A&& a,
			 const B& b,
			 const C& c)
  {
    const auto tmp=b[0]*c[0]-b[1]*c[1];
    
    a[1]=b[0]*c[1]+b[1]*c[0];
    a[0]=tmp;
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void complex_prodassign(A&& a,
			  const B& b)
  {
    safe_complex_prod(a,a,b);
  }
  
  //Minus version
  inline void safe_complex_prod_minus(complex a,const complex b,const complex c)
  {
    const double tmp=-(b[0]*c[0]-b[1]*c[1]);
    a[1]=-(b[0]*c[1]+b[1]*c[0]);
    a[0]=tmp;
  }
  
  //The product of a complex number by the conjugate of the second
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_complex_conj2_prod(A&& a,
			       const B& b,
			       const C& c)
  {
    const auto tmp=b[0]*c[0]+b[1]*c[1];
    
    a[1]=-b[0]*c[1]+b[1]*c[0];
    a[0]=tmp;
  }
  
  //Minus version
  inline void safe_complex_conj2_prod_minus(complex a,const complex b,const complex c)
  {
    const double tmp=-(b[0]*c[0]+b[1]*c[1]);
    a[1]=-(-b[0]*c[1]+b[1]*c[0]);
    a[0]=tmp;
  }
  
  //Swapped versions
  inline void safe_complex_conj1_prod(complex a,const complex b,const complex c)
  {safe_complex_conj2_prod(a,c,b);}
  inline void safe_complex_conj1_prod_minus(complex a,const complex b,const complex c)
  {safe_complex_conj2_prod_minus(a,c,b);}
  
  //complex prod i
  CUDA_HOST_AND_DEVICE inline void safe_complex_prod_i(complex a,const complex b)
  {
    const double tmp=b[0];
    a[0]=-b[1];
    a[1]=tmp;
  }
  CUDA_HOST_AND_DEVICE inline void assign_complex_prod_i(complex a)
  {safe_complex_prod_i(a,a);}
  
  //complex prod -i
  CUDA_HOST_AND_DEVICE inline void safe_complex_prod_minus_i(complex a,const complex b)
  {
    const double tmp=b[0];
    a[0]=b[1];
    a[1]=-tmp;
  }
  CUDA_HOST_AND_DEVICE inline void assign_complex_prod_minus_i(complex a)
  {safe_complex_prod_minus_i(a,a);}
  inline void complex_summ_the_prod_i(complex a,const complex b,const complex c)
  {
    a[1]+=b[0]*c[0]-b[1]*c[1];
    a[0]-=b[0]*c[1]+b[1]*c[0];
  }
  inline void complex_subt_the_prod_i(complex a,const complex b,const complex c)
  {
    a[1]-=b[0]*c[0]-b[1]*c[1];
    a[0]+=b[0]*c[1]+b[1]*c[0];
  }
  inline void complex_summ_the_conj2_prod_i(complex a,const complex b,const complex c)
  {
    a[1]+=+b[0]*c[0]+b[1]*c[1];
    a[0]-=-b[0]*c[1]+b[1]*c[0];
  }
  inline void complex_summ_the_conj1_prod_i(complex a,const complex b,const complex c)
  {complex_summ_the_conj2_prod(a,c,b);}
  inline void complex_summ_the_conj_conj_prod_i(complex a,const complex b,const complex c)
  {
    a[1]+=+b[0]*c[0]-b[1]*c[1];
    a[0]-=-b[0]*c[1]-b[1]*c[0];
  }
  inline void complex_subt_the_conj2_prod_i(complex a,const complex b,const complex c)
  {
    a[1]-=+b[0]*c[0]+b[1]*c[1];
    a[0]+=-b[0]*c[1]+b[1]*c[0];
  }
  inline void complex_subt_the_conj1_prod_i(complex a,const complex b,const complex c)
  {complex_subt_the_conj2_prod(a,c,b);}
  inline void complex_subt_the_conj_conj_prod_i(complex a,const complex b,const complex c)
  {
    a[1]-=+b[0]*c[0]-b[1]*c[1];
    a[0]+=-b[0]*c[1]-b[1]*c[0];
  }
  
  template <typename A>
  //squared norm
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  double complex_norm2(A&& c)
  {
    return c[0]*c[0]+c[1]*c[1];
  }
  
  //reciprocal of a complex
  CUDA_HOST_AND_DEVICE inline void complex_reciprocal(complex rec,const complex c)
  {
    const double module=c[0]*c[0]+c[1]*c[1];
    
    rec[0]=c[0]/module;
    rec[1]=-c[1]/module;
  }
  
  //squared root of a complex
  inline void complex_sqrt(complex res,const complex base)
  {
    double module=sqrt(base[0]*base[0]+base[1]*base[1]);
    const double cost=base[0]/module;
    const double sinth=sqrt(0.5*(1-cost));
    const double costh=sqrt(0.5*(1+cost));
    module=sqrt(module);
    if(base[1]>=0) res[0]=+module*costh;
    else           res[0]=-module*costh;
    res[1]=module*sinth;
  }
  
  //power of a complex
  CUDA_HOST_AND_DEVICE inline void complex_pow(complex res,const complex base,const double exp)
  {
    const double module=pow(base[0]*base[0]+base[1]*base[1],exp/2);
    const double anomaly=atan2(base[1],base[0])*exp;
    
    res[0]=module*cos(anomaly);
    res[1]=module*sin(anomaly);
  }
}
#endif
