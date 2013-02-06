#ifndef _FLOAT256_H
#define _FLOAT256_H
int float_256_is_equal(float_256 a,float_256 b);
int float_256_is_greater(float_256 a,double b);
int float_256_is_greater(float_256 a,float_256 b);
int float_256_is_smaller(float_256 a,double b);
int float_256_is_smaller(float_256 a,float_256 b);
void float_128_print(float_128 a);
void float_256_abs(float_256 a,float_256 b);
void float_256_copy(float_256 b,float_256 a);
void float_256_div(float_256 c,float_256 a,float_256 b);
void float_256_from_double(float_256 b,double a);
void float_256_pow_int(float_256 out,float_256 in,int d);
void float_256_pow_int_frac(float_256 out,float_256 ext_in,int n,int d);
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
#endif
