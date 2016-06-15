#ifndef _FLOAT_128_HPP
#define _FLOAT_128_HPP

#include "complex.hpp"

#ifndef EXTERN_FLOAT_128
 #define EXTERN_FLOAT_128 extern
#endif

#define NISSA_DEFAULT_USE_128_BIT_PRECISION 0

namespace nissa
{
  EXTERN_FLOAT_128 int use_128_bit_precision;
  
  //quadruple precision float
  typedef double float_128[2];
  typedef float_128 complex_128[2];
  typedef complex_128 color_128[NCOL];
  typedef color_128 spincolor_128[4];
  typedef complex_128 vir_complex_128[2];
  
  typedef vir_complex_128 vir_color_128[NCOL];
  typedef vir_color_128 vir_su3_128[NCOL];
  typedef vir_su3_128 vir_oct_su3_128[8];
  typedef vir_color_128 vir_spincolor_128[4];
  
  double double_from_float_128(float_128 b);
  void color_128_copy(color_128 a,color_128 b);
  void color_128_isubt(color_128 a,color_128 b,color_128 c);
  void color_128_isubtassign(color_128 a,color_128 b);
  void color_128_isumm(color_128 a,color_128 b,color_128 c);
  void color_128_isummassign(color_128 a,color_128 b);
  void color_128_subt(color_128 a,color_128 b,color_128 c);
  void color_128_subtassign(color_128 a,color_128 b);
  void color_128_summ(color_128 a,color_128 b,color_128 c);
  void color_128_summassign(color_128 a,color_128 b);
  void color_128_put_to_zero(color_128 a);
  void complex_128_isubt(complex_128 a,complex_128 b,complex_128 c);
  void complex_128_isumm(complex_128 a,complex_128 b,complex_128 c);
  void complex_128_subt(complex_128 a,complex_128 b,complex_128 c);
  void complex_128_summ(complex_128 a,complex_128 b,complex_128 c);
  void complex_128_summassign_64(complex_128 a,complex b);
  void complex_128_summassign(complex_128 a,complex_128 b);
  void complex_128_subtassign(complex_128 a,complex_128 b);
  void complex_summ_the_64_conj1_prod_128(complex_128 a,complex b,complex_128 c);
  void complex_summ_the_64_prod_128(complex_128 a,complex b,complex_128 c);
  void complex_subt_the_64_prod_128(complex_128 a,complex b,complex_128 c);
  void float_128_copy(float_128 b,float_128 a);
  void float_128_swap(float_128 b,float_128 a);
  int float_128_is_greater(float_128 a,float_128 b);
  int float_128_is_smaller(float_128 a,float_128 b);
  void float_128_abs(float_128 a,float_128 b);
  void float_128_put_to_zero(float_128 a);
  void float_128_from_64(float_128 b,double a);
  void float_128_print(float_128 a);
  void float_128_prod(float_128 c,float_128 a,float_128 b);
  void float_128_prodassign(float_128 out,float_128 in);
  void float_128_prod_64(float_128 c,float_128 a,double b);
  void float_128_64_prod_64(float_128 c,double a,double b);
  void float_128_pow_int_frac(float_128 out,float_128 in,int n,int d);
  void float_128_pow_int(float_128 out,float_128 in,int d);
  void float_128_div_64(float_128 div,float_128 a,double b);
  void float_128_div(float_128 div,float_128 a,float_128 b);
  void float_128_subt(float_128 c,float_128 a,float_128 b);
  void float_128_subt_from_64(float_128 c,double a,float_128 b);
  void float_128_subtassign(float_128 b,float_128 a);
  void float_128_summ(float_128 c,float_128 a,float_128 b);
  void float_128_64_summ_64(float_128 c,double a,double b);
  void float_128_summ_64(float_128 c,float_128 a,double b);
  void float_128_subt_the_prod(float_128 c,float_128 a,float_128 b);
  void float_128_summ_the_prod(float_128 c,float_128 a,float_128 b);
  void float_128_summ_the_64_prod(float_128 c,double a,double b);
  void float_128_summassign(float_128 b,float_128 a);
  void float_128_summassign_64(float_128 b,double a);
  void float_128_uminus(float_128 b,float_128 a);
  void float_64_prod_128(float_128 c,double a,float_128 b);
  void float_64_prod_complex_128(complex_128 a,double b,complex_128 c);
  void float_64_summ_the_iprod_complex_128(complex_128 a,double b,complex_128 c);
  void float_64_summ_the_prod_complex_128(complex_128 a,double b,complex_128 c);
  void float_subt_the_64_prod_128(float_128 c,double a,float_128 b);
  void float_summ_the_64_prod_128(float_128 c,double a,float_128 b);
  void unsafe_complex_64_prod_128(complex_128 a,complex b,complex_128 c);
  void unsafe_complex_64_conj1_prod_128(complex_128 a,complex b,complex_128 c);
}

#undef EXTERN_FLOAT_128

#endif
