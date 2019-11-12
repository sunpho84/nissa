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

namespace nissa
{
  typedef double complex[2];
  typedef complex quad_u1[NDIM];
  
  typedef float single_complex[2];
  
  typedef complex vir_complex[2];
  typedef single_complex vir_single_complex[2];
  
  //////////////////////////////////////////////////////////
  
  inline double real_part_of_complex_prod(const complex a,const complex b)
  {return a[0]*b[0]-a[1]*b[1];};
  CUDA_HOST_AND_DEVICE inline double real_part_of_complex_scalar_prod(const complex a,const complex b)
  {return a[0]*b[0]+a[1]*b[1];};
  //print
  inline void complex_print(const complex a)
  {printf("(%16.16lg,%16.16lg)\n",a[0],a[1]);}
  
  //Assign
  inline CUDA_HOST_AND_DEVICE void complex_copy(complex a,const complex b)
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
  CUDA_HOST_AND_DEVICE inline void complex_put_to_real(complex a,const double b)
  {
    a[0]=b;
    a[1]=0;
  }
  inline void complex_put_to_imag(complex a,const double b)
  {
    a[0]=0;
    a[1]=b;
  }
  CUDA_HOST_AND_DEVICE inline void complex_put_to_zero(complex a)
  {complex_put_to_real(a,0);}
  
  //Assign the conj
  CUDA_HOST_AND_DEVICE inline void complex_conj(complex a,const complex b)
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
  //The sum of two complex number
  CUDA_HOST_AND_DEVICE inline void complex_summ(complex a,const complex b,const complex c)
  {
    a[0]=b[0]+c[0];
    a[1]=b[1]+c[1];
  }
  CUDA_HOST_AND_DEVICE inline void complex_isumm(complex a,const complex b,const complex c)
  {
    a[0]=b[0]-c[1];
    a[1]=b[1]+c[0];
  }
  inline void complex_isummassign(complex a,const complex b)
  {complex_isumm(a,a,b);}
  CUDA_HOST_AND_DEVICE inline void complex_summ_conj2(complex a,const complex b,const complex c)
  {
    a[0]=b[0]+c[0];
    a[1]=b[1]-c[1];
  }
  inline void complex_summ_conj1(complex a,const complex b,const complex c)
  {complex_summ_conj2(a,c,b);}
  CUDA_HOST_AND_DEVICE inline void complex_subt(complex a,const complex b,const complex c)
  {
    a[0]=b[0]-c[0];
    a[1]=b[1]-c[1];
  }
  CUDA_HOST_AND_DEVICE inline void complex_isubt(complex a,const complex b,const complex c)
  {
    a[0]=b[0]+c[1];
    a[1]=b[1]-c[0];
  }
  inline void complex_isubtassign(complex a,const complex b)
  {complex_isubt(a,a,b);}
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
  CUDA_HOST_AND_DEVICE inline void complex_summassign(complex a,const complex b) {complex_summ(a,a,b);}
  inline void complex_subtassign(complex a,const complex b) {complex_subt(a,a,b);}
  
  //put to exp
  CUDA_HOST_AND_DEVICE inline void complex_iexp(complex out,const double arg)
  {sincos(arg,out+IM,out+RE);}
  
  //prod with real
  CUDA_HOST_AND_DEVICE inline void complex_prod_double(complex a,const complex b,const double c) {a[RE]=b[RE]*c;a[IM]=b[IM]*c;}
  inline void complex_prodassign_double(complex a,const double c) {complex_prod_double(a,a,c);}
  CUDA_HOST_AND_DEVICE inline void complex_prod_idouble(complex a,const complex b,const double c) {const double d=-b[IM]*c;a[IM]=b[RE]*c;a[RE]=d;}
  inline void complex_prodassign_idouble(complex a,const double b) {complex_prod_idouble(a,a,b);}
  
  //summ the prod with real
  CUDA_HOST_AND_DEVICE inline void complex_summ_the_prod_double(complex a,const complex b,const double c)
  {
    const double t=b[0]*c;
    a[1]+=b[1]*c;
    a[0]+=t;
  }
  inline void complex_subt_the_prod_double(complex a,const complex b,const double c)
  {
    const double t=b[0]*c;
    a[1]-=b[1]*c;
    a[0]-=t;
  }
  
  //summ the prod with imag
  CUDA_HOST_AND_DEVICE inline void complex_summ_the_prod_idouble(complex a,const complex b,double c)
  {
    const double t=b[1]*c;
    a[1]+=b[0]*c;
    a[0]-=t;
  }
  inline void complex_subt_the_prod_idouble(complex a,const complex b,double c)
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
  
  //Summ to the output the product of two complex number
  CUDA_HOST_AND_DEVICE inline void complex_summ_the_prod(complex a,const complex b,const complex c)
  {
    const double t=b[0]*c[0]-b[1]*c[1];
    a[1]+=b[0]*c[1]+b[1]*c[0];
    a[0]+=t;
  }
  CUDA_HOST_AND_DEVICE inline void single_complex_summ_the_prod(single_complex a,const single_complex b,const single_complex c)
  {
    const double t=b[0]*c[0]-b[1]*c[1];
    a[1]+=b[0]*c[1]+b[1]*c[0];
    a[0]+=t;
  }
  CUDA_HOST_AND_DEVICE inline void complex_subt_the_prod(complex a,const complex b,const complex c)
  {
    const double t=b[0]*c[0]-b[1]*c[1];
    a[1]-=b[0]*c[1]+b[1]*c[0];
    a[0]-=t;
  }
  CUDA_HOST_AND_DEVICE inline void complex_summ_the_conj2_prod(complex a,const complex b,const complex c)
  {
    const double t=+b[0]*c[0]+b[1]*c[1];
    a[1]+=-b[0]*c[1]+b[1]*c[0];
    a[0]+=t;
  }
  inline void single_complex_summ_the_conj2_prod(single_complex a,const single_complex b,const single_complex c)
  {
    const double t=+b[0]*c[0]+b[1]*c[1];
    a[1]+=-b[0]*c[1]+b[1]*c[0];
    a[0]+=t;
  }
  CUDA_HOST_AND_DEVICE inline void complex_summ_the_conj1_prod(complex a,const complex b,const complex c)
  {complex_summ_the_conj2_prod(a,c,b);}
  inline void single_complex_summ_the_conj1_prod(single_complex a,const single_complex b,const single_complex c)
  {single_complex_summ_the_conj2_prod(a,c,b);}
  CUDA_HOST_AND_DEVICE inline void complex_summ_the_conj_conj_prod(complex a,const complex b,const complex c)
  {
    const double t=+b[0]*c[0]-b[1]*c[1];
    a[1]+=-b[0]*c[1]-b[1]*c[0];
    a[0]+=t;
  }
  CUDA_HOST_AND_DEVICE inline void complex_subt_the_conj2_prod(complex a,const complex b,const complex c)
  {
    const double t=+b[0]*c[0]+b[1]*c[1];
    a[1]-=-b[0]*c[1]+b[1]*c[0];
    a[0]-=t;
  }
  CUDA_HOST_AND_DEVICE inline void single_complex_subt_the_conj2_prod(single_complex a,const single_complex b,const single_complex c)
  {
    const double t=+b[0]*c[0]+b[1]*c[1];
    a[1]-=-b[0]*c[1]+b[1]*c[0];
    a[0]-=t;
  }
  CUDA_HOST_AND_DEVICE inline void complex_subt_the_conj1_prod(complex a,const complex b,const complex c)
  {complex_subt_the_conj2_prod(a,c,b);}
  CUDA_HOST_AND_DEVICE inline void single_complex_subt_the_conj1_prod(single_complex a,const single_complex b,const single_complex c)
  {single_complex_subt_the_conj2_prod(a,c,b);}
  CUDA_HOST_AND_DEVICE inline void complex_subt_the_conj_conj_prod(complex a,const complex b,const complex c)
  {
    const double t=+b[0]*c[0]-b[1]*c[1];
    a[1]-=-b[0]*c[1]-b[1]*c[0];
    a[0]-=t;
  }
  
  //The product of two complex number
  CUDA_HOST_AND_DEVICE inline void unsafe_complex_prod(complex a,const complex b,const complex c)
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
  
  //The product of a complex number by the conjugate of the second
  CUDA_HOST_AND_DEVICE inline void unsafe_complex_conj2_prod(complex a,const complex b,const complex c)
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
  
  //Swapped order
  CUDA_HOST_AND_DEVICE inline void unsafe_complex_conj1_prod(complex a,const complex b,const complex c)
  {unsafe_complex_conj2_prod(a,c,b);}
  inline void unsafe_complex_conj1_prod_minus(complex a,const complex b,const complex c)
  {unsafe_complex_conj2_prod_minus(a,c,b);}
  
  //The product of the conjugate of two complex numbers
  CUDA_HOST_AND_DEVICE inline void unsafe_complex_conj_conj_prod(complex a,const complex b,const complex c)
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
  
  //The product of two complex number
  CUDA_HOST_AND_DEVICE inline void safe_complex_prod(complex a,const complex b,const complex c)
  {
    const double tmp=b[0]*c[0]-b[1]*c[1];
    a[1]=b[0]*c[1]+b[1]*c[0];
    a[0]=tmp;
  }
  inline void complex_prodassign(complex a,const complex b)
  {safe_complex_prod(a,a,b);}
  
  //Minus version
  inline void safe_complex_prod_minus(complex a,const complex b,const complex c)
  {
    const double tmp=-(b[0]*c[0]-b[1]*c[1]);
    a[1]=-(b[0]*c[1]+b[1]*c[0]);
    a[0]=tmp;
  }
  
  //The product of a complex number by the conjugate of the second
  inline void safe_complex_conj2_prod(complex a,const complex b,const complex c)
  {
    const double tmp=b[0]*c[0]+b[1]*c[1];
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
  //squared norm
  CUDA_HOST_AND_DEVICE inline double complex_norm2(const complex c)
  {return c[0]*c[0]+c[1]*c[1];}
  
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
