#ifndef _FLOAT256_HPP
#define _FLOAT256_HPP

#include "float_128.hpp"

namespace nissa
{
  //octpuple
  typedef double float_256[4];
  typedef double float_256_unr[5];
  struct float_256_class;
  
  int float_256_is_equal(float_256 a,float_256 b);
  int float_256_is_greater(float_256 a,double b);
  int float_256_is_greater(float_256 a,float_256 b);
  int float_256_is_smaller(float_256 a,double b);
  int float_256_is_smaller(float_256 a,float_256 b);
  //void float_128_print(float_128 a);
  void float_256_abs(float_256 a,float_256 b);
  void float_256_copy(float_256 b,float_256 a);
  void float_256_div(float_256 c,float_256 a,float_256 b);
  void float_256_from_double(float_256 b,double a);
  void float_256_pow_int(float_256 out,float_256 in,int d);
  void float_256_print(float_256 a);
  void float_256_prod(float_256 c,float_256 a,float_256 b);
  void float_256_prod_double(float_256 c,float_256 a,double b);
  void float_256_prodassign(float_256 out,float_256 in);
  void float_256_put_to_zero(float_256 a);
  void float_256_subt(float_256 c,float_256 a,float_256 b);
  void float_256_subt_the_prod(float_256 c,float_256 a,float_256 b);
  void float_256_subtassign(float_256 b,float_256 a);
  void float_256_summ(float_256 c,float_256 a,float_256 b);
  void float_256_summ_double(float_256 c,float_256 a,double b);
  void float_256_summ_ieee(float_256 c,float_256 a,float_256 b);
  void float_256_summ_the_prod(float_256 c,float_256 a,float_256 b);
  void float_256_summassign(float_256 b,float_256 a);
  void float_256_swap(float_256 b,float_256 a);
  void float_256_uminus(float_256 b,float_256 a);
  
  struct float_256_class
  {
    float_256 num;
    double &operator[](int i){return num[i];}
    double get_d(){return num[0];}
    
    float_256_class(double a){float_256_from_double(num,a);}
    float_256_class(){}
    float_256_class operator+(float_256_class a){float_256_class b;float_256_summ(b.num,this->num,a.num);return b;}
    float_256_class operator-(float_256_class a){float_256_class b;float_256_subt(b.num,this->num,a.num);return b;}
    float_256_class operator*(float_256_class a){float_256_class b;float_256_prod(b.num,this->num,a.num);return b;}
    float_256_class operator/(float_256_class a){float_256_class b;float_256_div(b.num,this->num,a.num);return b;}
    float_256_class operator+(double a){float_256_class b;float_256_summ_double(b.num,this->num,a);return b;}
    float_256_class operator*(double a){float_256_class b;float_256_prod_double(b.num,this->num,a);return b;}
    float_256_class operator-(double a){return operator+(-a);}
    float_256_class operator/(double a){return operator*(1/a);}
    float_256_class operator+=(float_256_class a){return (*this)=(*this)+a;}
    float_256_class operator-=(float_256_class a){return (*this)=(*this)-a;}
    float_256_class operator*=(float_256_class a){return (*this)=(*this)*a;}
    float_256_class operator/=(float_256_class a){return (*this)=(*this)/a;}
    float_256_class operator-(){float_256_class a;float_256_uminus(a.num,this->num);return a;}
    bool operator>(float_256_class a){return float_256_is_greater(this->num,a.num);}
    bool operator<(float_256_class a){return float_256_is_smaller(this->num,a.num);}
    bool operator>=(float_256_class a){return !(operator<(a));}
    bool operator<=(float_256_class a){return !(operator>(a));}
    bool operator==(float_256_class a){return float_256_is_equal(this->num,a.num);}
    bool operator!=(float_256_class a){return !(operator==(a));}
  };
  inline float_256_class operator+(double a,float_256_class b){return b+a;}
  inline float_256_class operator-(double a,float_256_class b){return (-b)+a;}
  inline float_256_class operator*(double a,float_256_class b){return b*a;}
  inline float_256_class operator/(double a,float_256_class b){return (float_256_class)a/b;}
  inline float_256_class abs(float_256_class a){float_256_class b;float_256_abs(b.num,a.num);return b;}
}

#endif
