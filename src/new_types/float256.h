#ifndef _FLOAT_256_H
#define _FLOAT_256_H
double float_256_most_sign_part(float_256 b);
int float_256_is_greater(float_256 a,float_256 b);
int float_256_is_smaller(float_256 a,float_256 b);
void float_128_from_256(float_128 a,float_256 b);
void float_128_prod_256(float_256 c,float_128 a,float_256 b);
void float_256_128_summ_128(float_256 c,float_128 a,float_128 b);
void float_256_abs(float_256 a,float_256 b);
void float_256_copy(float_256 b,float_256 a);
void float_256_div(float_256 div,float_256 a,float_256 b);
void float_256_div_128(float_256 div,float_256 a,float_128 b);
void float_256_from_128(float_256 b,float_128 a);
void float_256_from_double(float_256 b,double a);
void float_256_pow_int(float_256 out,float_256 in,int d);
void float_256_pow_int_frac(float_256 out,float_256 ext_in,int n,int d);
void float_256_print(float_256 a);
void float_256_prod(float_256 c,float_256 a,float_256 b);
void float_256_prod_128(float_256 c,float_256 a,float_128 b);
void float_256_prodassign(float_256 out,float_256 in);
void float_256_put_to_zero(float_256 a);
void float_256_subt(float_256 c,float_256 a,float_256 b);
void float_256_subt_from_128(float_256 c,float_128 a,float_256 b);
void float_256_subt_the_prod(float_256 c,float_256 a,float_256 b);
void float_256_subtassign(float_256 b,float_256 a);
void float_256_summ(float_256 c,float_256 a,float_256 b);
void float_256_summ_128(float_256 c,float_256 a,float_128 b);
void float_256_summ_the_128_prod(float_256 c,float_128 a,float_128 b);
void float_256_summ_the_prod(float_256 c,float_256 a,float_256 b);
void float_256_summassign(float_256 b,float_256 a);
void float_256_summassign_128(float_256 b,float_128 a);
void float_256_swap(float_256 b,float_256 a);
void float_256_uminus(float_256 b,float_256 a);
void float_subt_the_128_prod_256(float_256 c,float_128 a,float_256 b);
void float_summ_the_128_prod_256(float_256 c,float_128 a,float_256 b);
#endif
