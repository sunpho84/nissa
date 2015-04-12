#ifndef _COMPLEX_HPP
#define _COMPLEX_HPP

#include "new_types_definitions.hpp"

namespace nissa
{
  inline double real_part_of_complex_prod(complex a,complex b)
  {return a[0]*b[0]-a[1]*b[1];};
  inline double real_part_of_complex_scalar_prod(const complex a,const complex b)
  {return a[0]*b[0]+a[1]*b[1];};
  //print
  inline void complex_print(complex a)
  {printf("(%16.16lg,%16.16lg)\n",a[0],a[1]);}
  
  //Assign
  inline void complex_copy(complex a,complex b)
  {
    a[0]=b[0];
    a[1]=b[1];
  }
  inline void complex_copy_from_single_complex(complex a,single_complex b)
  {
    a[0]=b[0];
    a[1]=b[1];
  }
  inline void single_complex_copy_from_complex(single_complex a,complex b)
  {
    a[0]=b[0];
    a[1]=b[1];
  }
  inline void complex_put_to_real(complex a,double b)
  {
    a[0]=b;
    a[1]=0;
  }
  inline void complex_put_to_zero(complex a)
  {complex_put_to_real(a,0);}
  //Assign the conj
  inline void complex_conj(complex a,complex b)
  {
    a[0]=b[0];
    a[1]=-b[1];
  }
  //Assign minus the conj
  inline void complex_minus_conj(complex a,complex b)
  {
    a[0]=-b[0];
    a[1]=b[1];
  }
  //The sum of two complex number
  inline void complex_summ(complex a,complex b,complex c)
  {
    a[0]=b[0]+c[0];
    a[1]=b[1]+c[1];
  }
  inline void complex_isumm(complex a,complex b,complex c)
  {
    a[0]=b[0]-c[1];
    a[1]=b[1]+c[0];
  }
  inline void complex_isummassign(complex a,complex b)
  {complex_isumm(a,a,b);}
  inline void complex_summ_conj2(complex a,complex b,complex c)
  {
    a[0]=b[0]+c[0];
    a[1]=b[1]-c[1];
  }
  inline void complex_summ_conj1(complex a,complex b,complex c)
  {complex_summ_conj2(a,c,b);}
  inline void complex_subt(complex a,complex b,complex c)
  {
    a[0]=b[0]-c[0];
    a[1]=b[1]-c[1];
  }
  inline void complex_isubt(complex a,complex b,complex c)
  {
    a[0]=b[0]+c[1];
    a[1]=b[1]-c[0];
  }
  inline void complex_isubtassign(complex a,complex b)
  {complex_isubt(a,a,b);}
  inline void complex_subt_conj2(complex a,complex b,complex c)
  {
    a[0]=b[0]-c[0];
    a[1]=b[1]+c[1];
  }
  inline void complex_subt_conj1(complex a,complex b,complex c)
  {
    a[0]=b[0]-c[0];
    a[1]=-b[1]-c[1];
  }
  inline void complex_summassign(complex a,complex b)
  {
    a[0]+=b[0];
    a[1]+=b[1];
  }
  inline void complex_subtassign(complex a,complex b)
  {
    a[0]-=b[0];
    a[1]-=b[1];
  }
  
  //prod with real
  inline void complex_prod_double(complex a,complex b,double c)
  {
    a[0]=b[0]*c;
    a[1]=b[1]*c;
  }
  inline void complex_prodassign_double(complex a,double c)
  {complex_prod_double(a,a,c);}
  inline void unsafe_complex_prod_idouble(complex a,complex b,double c)
  {
    a[0]=-b[1]*c;
    a[1]= b[0]*c;
  }
  inline void safe_complex_prod_idouble(complex a,complex b,double c)
  {
    complex t;
    unsafe_complex_prod_idouble(t,b,c);
    complex_copy(a,t);
  }
  
  inline void complex_prodassign_idouble(complex a,double b)
  {
    double a0=a[0];
    a[0]=-a[1]*b;
    a[1]=   a0*b;
  }
  
  //summ the prod with real
  inline void complex_summ_the_prod_double(complex a,const complex b,const double c)
  {
    double t=b[0]*c;
    a[1]+=b[1]*c;
    a[0]+=t;
  }
  inline void complex_subt_the_prod_double(complex a,const complex b,const double c)
  {
    double t=b[0]*c;
    a[1]-=b[1]*c;
    a[0]-=t;
  }
  
  //summ the prod with imag
  inline void complex_summ_the_prod_idouble(complex a,complex b,double c)
  {
    double t=b[1]*c;
    a[1]+=b[0]*c;
    a[0]-=t;
  }
  inline void complex_subt_the_prod_idouble(complex a,complex b,double c)
  {
    double t=b[1]*c;
    a[1]-=b[0]*c;
    a[0]+=t;
  }
  
  //Summ to the output the product of two complex number
  inline void complex_summ_the_prod(complex a,complex b,complex c)
  {
    double t=b[0]*c[0]-b[1]*c[1];
    a[1]+=b[0]*c[1]+b[1]*c[0];
    a[0]+=t;
  }
  inline void single_complex_summ_the_prod(single_complex a,single_complex b,single_complex c)
  {
    double t=b[0]*c[0]-b[1]*c[1];
    a[1]+=b[0]*c[1]+b[1]*c[0];
    a[0]+=t;
  }
  inline void complex_subt_the_prod(complex a,complex b,complex c)
  {
    double t=b[0]*c[0]-b[1]*c[1];
    a[1]-=b[0]*c[1]+b[1]*c[0];
    a[0]-=t;
  }
  inline void complex_summ_the_conj2_prod(complex a,complex b,complex c)
  {
    double t=+b[0]*c[0]+b[1]*c[1];
    a[1]+=-b[0]*c[1]+b[1]*c[0];
    a[0]+=t;
  }
  inline void single_complex_summ_the_conj2_prod(single_complex a,single_complex b,single_complex c)
  {
    double t=+b[0]*c[0]+b[1]*c[1];
    a[1]+=-b[0]*c[1]+b[1]*c[0];
    a[0]+=t;
  }
  inline void complex_summ_the_conj1_prod(complex a,complex b,complex c)
  {complex_summ_the_conj2_prod(a,c,b);}
  inline void single_complex_summ_the_conj1_prod(single_complex a,single_complex b,single_complex c)
  {single_complex_summ_the_conj2_prod(a,c,b);}
  inline void complex_summ_the_conj_conj_prod(complex a,complex b,complex c)
  {
    double t=+b[0]*c[0]-b[1]*c[1];
    a[1]+=-b[0]*c[1]-b[1]*c[0];
    a[0]+=t;
  }
  inline void complex_subt_the_conj2_prod(complex a,complex b,complex c)
  {
    double t=+b[0]*c[0]+b[1]*c[1];
    a[1]-=-b[0]*c[1]+b[1]*c[0];
    a[0]-=t;
  }
  inline void single_complex_subt_the_conj2_prod(single_complex a,single_complex b,single_complex c)
  {
    double t=+b[0]*c[0]+b[1]*c[1];
    a[1]-=-b[0]*c[1]+b[1]*c[0];
    a[0]-=t;
  }
  inline void complex_subt_the_conj1_prod(complex a,complex b,complex c)
  {complex_subt_the_conj2_prod(a,c,b);}
  inline void single_complex_subt_the_conj1_prod(single_complex a,single_complex b,single_complex c)
  {single_complex_subt_the_conj2_prod(a,c,b);}
  inline void complex_subt_the_conj_conj_prod(complex a,complex b,complex c)
  {
    double t=+b[0]*c[0]-b[1]*c[1];
    a[1]-=-b[0]*c[1]-b[1]*c[0];
    a[0]-=t;
  }
  
  //The product of two complex number
  inline void unsafe_complex_prod(complex a,complex b,complex c)
  {
    a[0]=b[0]*c[0]-b[1]*c[1];
    a[1]=b[0]*c[1]+b[1]*c[0];
  }  
  inline void unsafe_single_single_complex_prod(single_complex a,single_complex b,single_complex c)
  {
    a[0]=b[0]*c[0]-b[1]*c[1];
    a[1]=b[0]*c[1]+b[1]*c[0];
  }  
  inline void unsafe_single_complex_prod(single_complex a,single_complex b,single_complex c)
  {
    a[0]=b[0]*c[0]-b[1]*c[1];
    a[1]=b[0]*c[1]+b[1]*c[0];
  }
  
  //Minus the product
  inline void unsafe_complex_prod_minus(complex a,complex b,complex c)
  {
    a[0]=-(b[0]*c[0]-b[1]*c[1]);
    a[1]=-(b[0]*c[1]+b[1]*c[0]);
  }
  
  //The product of a complex number by the conjugate of the second
  inline void unsafe_complex_conj2_prod(complex a,complex b,complex c)
  {
    a[0]=+b[0]*c[0]+b[1]*c[1];
    a[1]=-b[0]*c[1]+b[1]*c[0];
  }
  
  //Minus the former
  inline void unsafe_complex_conj2_prod_minus(complex a,complex b,complex c)
  {
    a[0]=-(b[0]*c[0]+b[1]*c[1]);
    a[1]=-(-b[0]*c[1]+b[1]*c[0]);
  }
  
  //Swapped order
  inline void unsafe_complex_conj1_prod(complex a,complex b,complex c)
  {unsafe_complex_conj2_prod(a,c,b);}
  inline void unsafe_complex_conj1_prod_minus(complex a,complex b,complex c)
  {unsafe_complex_conj2_prod_minus(a,c,b);}
  
  //The product of the conjugate of two complex numbers
  inline void unsafe_complex_conj_conj_prod(complex a,complex b,complex c)
  {
    a[0]=+b[0]*c[0]-b[1]*c[1];
    a[1]=-b[0]*c[1]-b[1]*c[0];
  }
  
  //Minus the former
  inline void unsafe_complex_conj_conj_prod_minus(complex a,complex b,complex c)
  {
    a[0]=-b[0]*c[0]+b[1]*c[1];
    a[1]=+b[0]*c[1]+b[1]*c[0];
  }
  
  //The product of two complex number
  inline void safe_complex_prod(complex a,complex b,complex c)
  {
    double tmp=b[0]*c[0]-b[1]*c[1];
    a[1]=b[0]*c[1]+b[1]*c[0];
    a[0]=tmp;
  }
  
  //Minus version
  inline void safe_complex_prod_minus(complex a,complex b,complex c)
  {
    double tmp=-(b[0]*c[0]-b[1]*c[1]);
    a[1]=-(b[0]*c[1]+b[1]*c[0]);
    a[0]=tmp;
  }
  
  //The product of a complex number by the conjugate of the second
  inline void safe_complex_conj2_prod(complex a,complex b,complex c)
  {
    double tmp=b[0]*c[0]+b[1]*c[1];
    a[1]=-b[0]*c[1]+b[1]*c[0];
    a[0]=tmp;
  }
  
  //Minus version
  inline void safe_complex_conj2_prod_minus(complex a,complex b,complex c)
  {
    double tmp=-(b[0]*c[0]+b[1]*c[1]);
    a[1]=-(-b[0]*c[1]+b[1]*c[0]);
    a[0]=tmp;
  }
  
  //Swapped versions
  inline void safe_complex_conj1_prod(complex a,complex b,complex c)
  {safe_complex_conj2_prod(a,c,b);}
  inline void safe_complex_conj1_prod_minus(complex a,complex b,complex c)
  {safe_complex_conj2_prod_minus(a,c,b);}
  
  //complex prod i
  inline void safe_complex_prod_i(complex a,complex b)
  {
    double tmp=b[0];
    a[0]=-b[1];
    a[1]=tmp;
  }
  inline void assign_complex_prod_i(complex a)
  {safe_complex_prod_i(a,a);}
  
  //complex prod -i
  inline void safe_complex_prod_minus_i(complex a,complex b)
  {
    double tmp=b[0];
    a[0]=b[1];
    a[1]=-tmp;
  }
  inline void assign_complex_prod_minus_i(complex a)
  {safe_complex_prod_minus_i(a,a);}
  inline void complex_summ_the_prod_i(complex a,complex b,complex c)
  {
    a[1]+=b[0]*c[0]-b[1]*c[1];
    a[0]-=b[0]*c[1]+b[1]*c[0];
  }
  inline void complex_subt_the_prod_i(complex a,complex b,complex c)
  {
    a[1]-=b[0]*c[0]-b[1]*c[1];
    a[0]+=b[0]*c[1]+b[1]*c[0];
  }
  inline void complex_summ_the_conj2_prod_i(complex a,complex b,complex c)
  {
    a[1]+=+b[0]*c[0]+b[1]*c[1];
    a[0]-=-b[0]*c[1]+b[1]*c[0];
  }
  inline void complex_summ_the_conj1_prod_i(complex a,complex b,complex c)
  {complex_summ_the_conj2_prod(a,c,b);}
  inline void complex_summ_the_conj_conj_prod_i(complex a,complex b,complex c)
  {
    a[1]+=+b[0]*c[0]-b[1]*c[1];
    a[0]-=-b[0]*c[1]-b[1]*c[0];
  }
  inline void complex_subt_the_conj2_prod_i(complex a,complex b,complex c)
  {
    a[1]-=+b[0]*c[0]+b[1]*c[1];
    a[0]+=-b[0]*c[1]+b[1]*c[0];
  }
  inline void complex_subt_the_conj1_prod_i(complex a,complex b,complex c)
  {complex_subt_the_conj2_prod(a,c,b);}
  inline void complex_subt_the_conj_conj_prod_i(complex a,complex b,complex c)
  {
    a[1]-=+b[0]*c[0]-b[1]*c[1];
    a[0]+=-b[0]*c[1]-b[1]*c[0];
  }
  //squared norm
  inline double squared_complex_norm(complex c)
  {return c[0]*c[0]+c[1]*c[1];}
  
  //reciprocal of a complex
  inline void complex_reciprocal(complex rec,complex c)
  {
    double module=c[0]*c[0]+c[1]*c[1];
    
    rec[0]=c[0]/module;
    rec[1]=-c[1]/module;
  }
  
  //squared root of a complex
  inline void complex_sqrt(complex res,complex base)
  {
    double module=sqrt(base[0]*base[0]+base[1]*base[1]);
    double cost=base[0]/module;
    double sinth=sqrt(0.5*(1-cost));
    double costh=sqrt(0.5*(1+cost));
    module=sqrt(module);
    if(base[1]>=0) res[0]=+module*costh;
    else           res[0]=-module*costh;
    res[1]=module*sinth;
  }
  
  //power of a complex
  inline void complex_pow(complex res,complex base,double exp)
  {
    double module=pow(base[0]*base[0]+base[1]*base[1],exp/2);
    double anomaly=atan2(base[1],base[0])*exp;
    
    res[0]=module*cos(anomaly);
    res[1]=module*sin(anomaly);
  }
}
#endif
  
