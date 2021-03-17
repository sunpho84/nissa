#ifndef _FLOAT_128_HPP
#define _FLOAT_128_HPP

#ifndef EXTERN_FLOAT_128
 #define EXTERN_FLOAT_128 extern
#endif

#include <array>

#include "base/debug.hpp"
#include "new_types/dirac.hpp"
#include "new_types/su3.hpp"

#include "complex.hpp"

#define NISSA_DEFAULT_USE_128_BIT_PRECISION 0

#if defined(__ICC)
#pragma optimize("", off)
#endif

namespace nissa
{
  EXTERN_FLOAT_128 int use_128_bit_precision;
  
  //quadruple precision float
  struct float_128
  {
    double data[2];
    
    CUDA_HOST_AND_DEVICE double& operator[](const int i)
    {
      return data[i];
    }
    
    CUDA_HOST_AND_DEVICE const double& operator[](const int i) const
    {
      return data[i];
    }
  };
  
  typedef float_128 complex_128[2];
  typedef complex_128 color_128[NCOL];
  typedef color_128 halfspincolor_128[NDIRAC/2];
  typedef color_128 spincolor_128[NDIRAC];
  typedef complex_128 vir_complex_128[2];
  
  typedef vir_complex_128 vir_color_128[NCOL];
  typedef vir_color_128 vir_su3_128[NCOL];
  typedef vir_su3_128 vir_oct_su3_128[8];
  typedef vir_color_128 vir_spincolor_128[NDIRAC];
  
  inline void float_128_copy(float_128& b,const float_128 a)
  {
    b[0]=a[0];
    b[1]=a[1];
  }
  
  inline void float_128_swap(float_128& b,float_128& a)
  {
    float_128 c;
    float_128_copy(c,b);
    float_128_copy(b,a);
    float_128_copy(a,c);
  }
  
  CUDA_HOST_AND_DEVICE inline void float_128_uminus(float_128& b,const float_128& a)
  {
    b[0]=-a[0];
    b[1]=-a[1];
  }
  
  CUDA_HOST_AND_DEVICE inline void float_128_from_64(float_128& b,const double& a)
  {
    b[0]=a;
    b[1]=0;
  }
  
  CUDA_HOST_AND_DEVICE inline void float_128_put_to_zero(float_128& a)
  {a[0]=a[1]=0;}
  
  CUDA_HOST_AND_DEVICE inline double double_from_float_128(float_128 b)
  {return b[0]+b[1];}
  
  //128 summ 128
  CUDA_HOST_AND_DEVICE inline void float_128_summ(float_128& c,const float_128& a,const float_128& b)
  {
#ifndef fake_128
    double t1=a[0]+b[0];
    double e=t1-a[0];
    double t2=((b[0]-e)+(a[0]-(t1-e)))+a[1]+b[1];
    
    c[0]=t1+t2;
    c[1]=t2-(c[0]-t1);
#else
    c[0]=a[0]+b[0];
    c[1]=0;
#endif
  }
  CUDA_HOST_AND_DEVICE inline void float_128_summassign(float_128& b,const float_128& a)
  {float_128_summ(b,b,a);}
  CUDA_HOST_AND_DEVICE inline void float_128_subt(float_128& c,const float_128& a,const float_128& b)
  {
    float_128 d;
    float_128_uminus(d,b);
    float_128_summ(c,d,a);
  }
  CUDA_HOST_AND_DEVICE inline void float_128_subtassign(float_128& b,const float_128& a)
  {float_128_subt(b,b,a);}
  
  //128 summ 64
  CUDA_HOST_AND_DEVICE inline void float_128_summ_64(float_128& c,const float_128& a,const double& b)
  {
#ifndef fake_128
    double t1=a[0]+b;
    double e=t1-a[0];
    double t2=((b-e)+(a[0]-(t1-e)))+a[1];
    
    c[0]=t1+t2;
    c[1]=t2-(c[0]-t1);
#else
    c[0]=a[0]+b;
    c[1]=0;
#endif
  }
  
  //64 summ 64
  inline void float_128_64_summ_64(float_128& c,const double& a,const double& b)
  {
#ifndef fake_128
    double t1=a+b;
    double e=t1-a;
    double t2=((b-e)+(a-(t1-e)));
    
    c[0]=t1+t2;
    c[1]=t2-(c[0]-t1);
#else
    c[0]=a+b;
    c[1]=0;
#endif
  }
  
  CUDA_HOST_AND_DEVICE inline void float_128_summassign_64(float_128& b,const double& a)
  {
    float_128_summ_64(b,b,a);
  }
  
  CUDA_HOST_AND_DEVICE inline void float_128_subt_from_64(float_128& c,const double& a,const float_128& b)
  {
    float_128 d;
    float_128_uminus(d,b);
    float_128_summ_64(c,d,a);
  }
  
  inline void float_128_subtassign_64(float_128& b,const double& a)
  {
    float_128_summassign_64(b,-a);
  }
  
  //128 prod 128
  inline void float_128_prod(float_128& c,const float_128& a,const float_128& b)
  {
#ifndef fake_128
    const double split=134217729;
    
    double cona=a[0]*split;
    double conb=b[0]*split;
    
    double a1=cona-(cona-a[0]);
    double b1=conb-(conb-b[0]);
    double a2=a[0]-a1;
    double b2=b[0]-b1;
    
    double c11=a[0]*b[0];
    double c21=a2*b2+(a2*b1+(a1*b2+(a1*b1-c11)));
    
    double c2=a[0]*b[1]+a[1]*b[0];
    
    double t1=c11+c2;
    double e=t1-c11;
    double t2=a[1]*b[1]+((c2-e)+(c11-(t1-e)))+c21;
    
    c[0]=t1+t2;
    c[1]=t2-(c[0]-t1);
#else
    c[0]=a[0]*b[0];
    c[1]=0;
#endif
  }
  inline void float_128_summ_the_prod(float_128& c,const float_128& a,const float_128& b)
  {
    float_128 d;
    float_128_prod(d,a,b);
    float_128_summassign(c,d);
  }
  inline void float_128_subt_the_prod(float_128& c,const float_128& a,const float_128& b)
  {
    float_128 d;
    float_128_prod(d,a,b);
    float_128_subtassign(c,d);
  }
  
  inline void float_128_summ_the_64_prod(float_128& c,const double& a,const double& b)
  {
    float_128_summassign_64(c,a*b);
  }
  
  inline void float_128_subt_the_64_prod(float_128& c,const double& a,const double& b)
  {
    float_128_subtassign_64(c,a*b);
  }
  
  //128 prod 64
  CUDA_HOST_AND_DEVICE inline void float_128_prod_64(float_128& c,const float_128& a,const double& b)
  {
#ifndef fake_128
    const double split=134217729;
    
    double cona=a[0]*split;
    double conb=b*split;
    
    double a1=cona-(cona-a[0]);
    double b1=conb-(conb-b);
    double a2=a[0]-a1;
    double b2=b-b1;
    
    double c11=a[0]*b;
    double c21=a2*b2+(a2*b1+(a1*b2+(a1*b1-c11)));
    
    double c2=a[1]*b;
    
    double t1=c11+c2;
    double e=t1-c11;
    double t2=((c2-e)+(c11-(t1-e)))+c21;
    
    c[0]=t1+t2;
    c[1]=t2-(c[0]-t1);
#else
    c[0]=a[0]*b;
    c[1]=0;
#endif
  }
  CUDA_HOST_AND_DEVICE inline void float_64_prod_128(float_128& c,const double& a,const float_128& b)
  {float_128_prod_64(c,b,a);}
  CUDA_HOST_AND_DEVICE inline void float_summ_the_64_prod_128(float_128& c,const double& a,const float_128& b)
  {
    float_128 d;
    float_64_prod_128(d,a,b);
    float_128_summassign(c,d);
  }
  CUDA_HOST_AND_DEVICE inline void float_subt_the_64_prod_128(float_128& c,const double& a,const float_128& b)
  {
    float_128 d;
    float_64_prod_128(d,a,b);
    float_128_subtassign(c,d);
  }
  inline void float_128_64_prod_64(float_128& c,const double& a,const double& b)
  {
#ifndef fake_128
    const double split=134217729;
    
    double cona=a*split;
    double conb=b*split;
    
    double a1=cona-(cona-a);
    double b1=conb-(conb-b);
    double a2=a-a1;
    double b2=b-b1;
    
    double c11=a*b;
    double c21=a2*b2+(a2*b1+(a1*b2+(a1*b1-c11)));
    
    c[0]=c11+c21;
    c[1]=c21-(c[0]-c11);
#else
    c[0]=a*b;
    c[1]=0;
#endif
  }
  inline void float_128_prodassign(float_128& out,const float_128& in)
  {float_128_prod(out,out,in);}
  
  //divide two float_128
  inline void float_128_div(float_128& div,const float_128& a,const float_128& b)
  {
    //compute dividing factor
    double c=1/(b[0]+b[1]);
    //compute approx div
    float_128 div1;
    float_128_prod_64(div1,a,c);
    //compute first remaining
    float_128 rem;
    float_128_prod(rem,div1,b);
    float_128_subt(rem,a,rem);
    //compute the second piece
    float_128 div2;
    float_128_prod_64(div2,rem,c);
    //summ the two pieces
    float_128 div12;
    float_128_summ(div12,div1,div2);
    //second remaining
    float_128_prod(rem,div12,b);
    float_128_subt(rem,a,rem);
    //compute the third piece
    float_128 div3;
    float_128_prod_64(div3,rem,c);
    //summ the two pieces
    float_128_summ(div,div12,div3);
  }
  
  //divide a float_128 by a double
  inline void float_128_div_64(float_128& div,const float_128& a,const double& b)
  {
    //compute dividing factor
    double c=1/b;
    //compute approx div
    float_128 div1;
    float_128_prod_64(div1,a,c);
    //compute first remaining
    float_128 rem;
    float_128_prod_64(rem,div1,b);
    float_128_subt(rem,a,rem);
    //compute the second piece
    float_128 div2;
    float_128_prod_64(div2,rem,c);
    //summ the two pieces
    float_128 div12;
    float_128_summ(div12,div1,div2);
    //second remaining
    float_128_prod_64(rem,div12,b);
    float_128_subt(rem,a,rem);
    //compute the third piece
    float_128 div3;
    float_128_prod_64(div3,rem,c);
    //summ the two pieces
    float_128_summ(div,div12,div3);
  }
  
  //integer power
  inline void float_128_pow_int(float_128& out,const float_128& in,const int& d)
  {
    //negative or null case
    if(d<=0)
      {
	float_128_from_64(out,1);
	//if negative
	if(d<0)
	  {
	    //compute inv and assign to out
	    float_128 inv;
	    float_128_div(inv,out,in);
	    float_128_copy(out,inv);
	    //multiply the remaining numer of times
	    for(int i=2;i<=-d;i++) float_128_prodassign(out,inv);
	  }
      }
    //positive case
    else
      {
	float_128_copy(out,in);
	for(int i=2;i<=d;i++) float_128_prodassign(out,in);
      }
  }
  
  //frac power
  inline void float_128_pow_int_frac(float_128& out,const float_128& in,const int& n,const int& d)
  {
    //compute by solving out^d=in^n=ref
    float_128 ref;
    float_128_pow_int(ref,in,n);
    
    //let's start from a reasonable approx
    double r=(double)n/d;
    double sto=pow(in[0],r-1);
    float_128_64_summ_64(out,sto*in[0],sto*r*in[1]);
    //another possible approach
    //float_128_from_64(out,1+(in[0]-1)*r);
    
    //(out+err)^d=in^n -> err=out*rel_err, rel_err=(ref/out^d-1)/d
    float_128 rel_err;
    do
      {
	//compute out^d
	float_128 outd;
	float_128_pow_int(outd,out,d);
	
	//divide ref by out^d and subtract 1
	float_128_div(rel_err,ref,outd);
	float_128_summassign_64(rel_err,-1);
	float_128_prod_64(rel_err,rel_err,1.0/d);
	
	//printf("rel err: ");
	//float_128_print(rel_err);
	//printf("\n");
	
	//total err
	float_128 err;
	float_128_prod(err,rel_err,out);
	
	float_128_summassign(out,err);
      }
    while(fabs(rel_err[0])>3.e-32);
  }
  
  //a>b?
  inline int float_128_is_greater(const float_128& a,const float_128& b)
  {
    if(a[0]>b[0]) return true;
    if(a[0]<b[0]) return false;
    
    if(a[1]>b[1]) return true;
    
    return false;
  }
  
  //a<b?
  inline int float_128_is_smaller(const float_128& a,const float_128& b)
  {
    if(a[0]<b[0]) return true;
    if(a[0]>b[0]) return false;
    
    if(a[1]<b[1]) return true;
    
    return false;
  }
  
  inline void float_128_abs(float_128& a,const float_128& b)
  {
    if(b[0]>0||(b[0]==0&&b[1]>0)) float_128_copy(a,b);
    else float_128_uminus(a,b);
  }
  
  //////////////////////////////////////////////////////
  
  CUDA_HOST_AND_DEVICE inline void complex_128_put_to_zero(complex_128& a)
  {for(int ri=0;ri<2;ri++) float_128_put_to_zero(a[ri]);}
  
  //c128 summ c128
  CUDA_HOST_AND_DEVICE inline void complex_128_summ(complex_128& a,const complex_128& b,const complex_128& c)
  {for(int ri=0;ri<2;ri++) float_128_summ(a[ri],b[ri],c[ri]);}
  CUDA_HOST_AND_DEVICE inline void complex_128_summassign(complex_128& a,const complex_128& b)
  {complex_128_summ(a,a,b);}
  inline void complex_128_summassign_64(complex_128& a,const complex& b)
  {for(int ri=0;ri<2;ri++) float_128_summassign_64(a[ri],b[ri]);}
  CUDA_HOST_AND_DEVICE inline void complex_128_subt(complex_128& a,const complex_128& b,const complex_128& c)
  {for(int ri=0;ri<2;ri++) float_128_subt(a[ri],b[ri],c[ri]);}
  CUDA_HOST_AND_DEVICE inline void complex_128_subtassign(complex_128& a,const complex_128& b)
  {complex_128_subt(a,a,b);}
  
  //c128 isumm c128
  CUDA_HOST_AND_DEVICE inline void complex_128_isumm(complex_128& a,const complex_128& b,const complex_128& c)
  {
    float_128 d={c[0][0],c[0][1]};
    float_128_subt(a[0],b[0],c[1]);
    float_128_summ(a[1],b[1],d);
  }
  CUDA_HOST_AND_DEVICE inline void complex_128_isubt(complex_128& a,const complex_128& b,const complex_128& c)
  {
    float_128 d={c[0][0],c[0][1]};
    float_128_summ(a[0],b[0],c[1]);
    float_128_subt(a[1],b[1],d);
  }
  
  //64 prod c128
  CUDA_HOST_AND_DEVICE inline void float_64_prod_complex_128(complex_128& a,const double& b,const complex_128& c)
  {
    float_64_prod_128(a[0],b,c[0]);
    float_64_prod_128(a[1],b,c[1]);
  }
  CUDA_HOST_AND_DEVICE inline void float_64_summ_the_prod_complex_128(complex_128& a,const double& b,const complex_128& c)
  {
    float_summ_the_64_prod_128(a[0],b,c[0]);
    float_summ_the_64_prod_128(a[1],b,c[1]);
  }
  
  //64 iprod c128
  CUDA_HOST_AND_DEVICE inline void float_64_summ_the_iprod_complex_128(complex_128& a,const double& b,const complex_128& c)
  {
    float_128 d={c[0][0],c[0][1]};
    float_subt_the_64_prod_128(a[0],b,c[1]);
    float_summ_the_64_prod_128(a[1],b,d);
  }
  
  //c64 prod c128
  CUDA_HOST_AND_DEVICE inline void unsafe_complex_64_prod_128(complex_128& a,const complex& b,const complex_128& c)
  {
    //real part
    float_64_prod_128(a[0],b[0],c[0]);
    float_subt_the_64_prod_128(a[0],b[1],c[1]);
    //imag part
    float_64_prod_128(a[1],b[0],c[1]);
    float_summ_the_64_prod_128(a[1],b[1],c[0]);
  }
  CUDA_HOST_AND_DEVICE inline void unsafe_complex_128_prod_64(complex_128& a,const complex_128& b,const complex& c)
  {unsafe_complex_64_prod_128(a,c,b);}
  CUDA_HOST_AND_DEVICE inline void complex_summ_the_64_prod_128(complex_128& a,const complex& b,const complex_128& c)
  {
    complex_128 d;
    unsafe_complex_64_prod_128(d,b,c);
    complex_128_summassign(a,d);
  }
  CUDA_HOST_AND_DEVICE inline void complex_subt_the_64_prod_128(complex_128& a,const complex& b,const complex_128& c)
  {
    complex_128 d;
    unsafe_complex_64_prod_128(d,b,c);
    complex_128_subtassign(a,d);
  }
  
  //c64~ prod c128
  CUDA_HOST_AND_DEVICE inline void unsafe_complex_64_conj1_prod_128(complex_128& a,const complex& b,const complex_128& c)
  {
    complex d;
    complex_conj(d,b);
    unsafe_complex_64_prod_128(a,d,c);
  }
  CUDA_HOST_AND_DEVICE inline void complex_summ_the_64_conj1_prod_128(complex_128& a,const complex& b,const complex_128& c)
  {
    complex d;
    complex_conj(d,b);
    complex_summ_the_64_prod_128(a,d,c);
  }
  
  inline void color_128_put_to_zero(color_128& a)
  {memset(a,0,sizeof(color_128));}
  CUDA_HOST_AND_DEVICE inline void color_128_copy(color_128& a,const color_128& b)
  {memcpy(a,b,sizeof(color_128));}
  
  CUDA_HOST_AND_DEVICE inline void color_128_summ(color_128& a,const color_128& b,const color_128& c)
  {for(int ic=0;ic<3;ic++) complex_128_summ(a[ic],b[ic],c[ic]);}
  CUDA_HOST_AND_DEVICE inline void color_128_summassign(color_128& a,const color_128& b)
  {color_128_summ(a,a,b);}
  
  CUDA_HOST_AND_DEVICE inline void color_128_isumm(color_128& a,const color_128& b,const color_128& c)
  {for(int ic=0;ic<3;ic++) complex_128_isumm(a[ic],b[ic],c[ic]);}
  CUDA_HOST_AND_DEVICE inline void color_128_isummassign(color_128& a,const color_128& b)
  {color_128_isumm(a,a,b);}
  
  CUDA_HOST_AND_DEVICE inline void color_128_subt(color_128& a,const color_128& b,const color_128& c)
  {for(int ic=0;ic<3;ic++) complex_128_subt(a[ic],b[ic],c[ic]);}
  CUDA_HOST_AND_DEVICE inline void color_128_subtassign(color_128& a,const color_128& b)
  {color_128_subt(a,a,b);}
  
  CUDA_HOST_AND_DEVICE inline void color_128_isubt(color_128& a,const color_128& b,const color_128& c)
  {for(int ic=0;ic<3;ic++) complex_128_isubt(a[ic],b[ic],c[ic]);}
  CUDA_HOST_AND_DEVICE inline void color_128_isubtassign(color_128& a,const color_128& b)
  {color_128_isubt(a,a,b);}
  
  CUDA_HOST_AND_DEVICE inline void unsafe_color_128_prod_complex_64(color_128& out,const color_128& in,const complex& factor)
  {for(size_t i=0;i<NCOL;i++) unsafe_complex_128_prod_64(((complex_128*)out)[i],((complex_128*)in)[i],factor);}
  
  CUDA_HOST_AND_DEVICE inline void unsafe_su3_prod_color_128(color_128& a,const su3& b,const color_128& c)
  {
    for(int c1=0;c1<NCOL;c1++)
      {
	unsafe_complex_64_prod_128(a[c1],b[c1][0],c[0]);
	for(int c2=1;c2<NCOL;c2++) complex_summ_the_64_prod_128(a[c1],b[c1][c2],c[c2]);
      }
  }
  
  CUDA_HOST_AND_DEVICE inline void unsafe_su3_dag_prod_color_128(color_128& a,const su3& b,const color_128& c)
  {
    for(int c1=0;c1<NCOL;c1++)
      {
	unsafe_complex_64_conj1_prod_128(a[c1],b[0][c1],c[0]);
	for(int c2=1;c2<NCOL;c2++) complex_summ_the_64_conj1_prod_128(a[c1],b[c2][c1],c[c2]);
      }
  }
  
  CUDA_HOST_AND_DEVICE inline void su3_dag_summ_the_prod_color_128(color_128& a,const su3& b,const color_128& c)
  {
    for(int c1=0;c1<NCOL;c1++)
      for(int c2=0;c2<NCOL;c2++)
	complex_summ_the_64_conj1_prod_128(a[c1],b[c2][c1],c[c2]);
  }
  
  CUDA_HOST_AND_DEVICE inline void su3_subt_the_prod_color_128(color_128& a,const su3& b,const color_128& c)
  {
    for(int c1=0;c1<NCOL;c1++)
      for(int c2=0;c2<NCOL;c2++)
	complex_subt_the_64_prod_128(a[c1],b[c1][c2],c[c2]);
  }
  
  CUDA_HOST_AND_DEVICE inline void su3_summ_the_prod_color_128(color_128& a,const su3& b,const color_128& c)
  {
    for(int c1=0;c1<NCOL;c1++)
      for(int c2=0;c2<NCOL;c2++)
	complex_summ_the_64_prod_128(a[c1],b[c1][c2],c[c2]);
  }
  
  CUDA_HOST_AND_DEVICE inline void unsafe_halfspincolor_halfspincolor_times_halfspincolor_128(halfspincolor_128& a,const halfspincolor_halfspincolor& b,const halfspincolor_128& c)
  {
    for(int id_out=0;id_out<NDIRAC/2;id_out++)
      for(int ic_out=0;ic_out<NCOL;ic_out++)
	{
	  complex_128_put_to_zero(a[id_out][ic_out]);
	  for(int id_in=0;id_in<NDIRAC/2;id_in++)
	    for(int ic_in=0;ic_in<NCOL;ic_in++)
	      complex_summ_the_64_prod_128(a[id_out][ic_out],b[id_out][ic_out][id_in][ic_in],c[id_in][ic_in]);
	}
  }
  
  CUDA_HOST_AND_DEVICE inline void unsafe_halfspincolor_halfspincolor_dag_times_halfspincolor_128(halfspincolor_128& a,const halfspincolor_halfspincolor& b,const halfspincolor_128& c)
  {
    for(int id_out=0;id_out<NDIRAC/2;id_out++)
      for(int ic_out=0;ic_out<NCOL;ic_out++)
	{
	  complex_128_put_to_zero(a[id_out][ic_out]);
	  for(int id_in=0;id_in<NDIRAC/2;id_in++)
	    for(int ic_in=0;ic_in<NCOL;ic_in++)
	      complex_summ_the_64_conj1_prod_128(a[id_out][ic_out],b[id_in][ic_in][id_out][ic_out],c[id_in][ic_in]);
	}
  }
  
  template <>
  inline float_128 summ(const float_128& a,const float_128& b)
  {
    float_128 c;
    float_128_summ(c,a,b);
    
    return c;
  }
  
}

#undef EXTERN_FLOAT_128

#if defined(__ICC)
#pragma optimize("", on)
#endif

#endif
