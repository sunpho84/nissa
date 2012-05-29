#ifndef _COMPLEX_H
#define _COMPLEX_H

#include "new_types_definitions.h"

double squared_complex_norm(complex c);
void assign_complex_prod_i(complex a);
void assign_complex_prod_minus_i(complex a);
void complex_conj(complex a,complex b);
void complex_copy(complex a,complex b);
void complex_isubt(complex a,complex b,complex c);
void complex_isubtassign(complex a,complex b);
void complex_isumm(complex a,complex b,complex c);
void complex_isummassign(complex a,complex b);
void complex_minus_conj(complex a,complex b);
void complex_pow(complex res,complex base,double exp);
void complex_prod_double(complex a,complex b,double c);
void complex_prodassign_double(complex a,double c);
void complex_prodassign_idouble(complex a,double c);
void complex_reciprocal(complex rec,complex c);
void complex_sqrt(complex res,complex base);
void complex_subt(complex a,complex b,complex c);
void complex_subt_conj1(complex a,complex b,complex c);
void complex_subt_conj2(complex a,complex b,complex c);
void complex_subt_the_conj1_prod(complex a,complex b,complex c);
void complex_subt_the_conj1_prod_i(complex a,complex b,complex c);
void complex_subt_the_conj2_prod(complex a,complex b,complex c);
void complex_subt_the_conj2_prod_i(complex a,complex b,complex c);
void complex_subt_the_conj_conj_prod(complex a,complex b,complex c);
void complex_subt_the_conj_conj_prod_i(complex a,complex b,complex c);
void complex_subt_the_prod(complex a,complex b,complex c);
void complex_subt_the_prod_double(complex a,complex b,double c);
void complex_subt_the_prod_i(complex a,complex b,complex c);
void complex_subt_the_prod_idouble(complex a,complex b,double c);
void complex_subtassign(complex a,complex b);
void complex_summ(complex a,complex b,complex c);
void complex_summ_conj1(complex a,complex b,complex c);
void complex_summ_conj2(complex a,complex b,complex c);
void complex_summ_the_conj1_prod(complex a,complex b,complex c);
void complex_summ_the_conj1_prod_i(complex a,complex b,complex c);
void complex_summ_the_conj2_prod(complex a,complex b,complex c);
void complex_summ_the_conj2_prod_i(complex a,complex b,complex c);
void complex_summ_the_conj_conj_prod(complex a,complex b,complex c);
void complex_summ_the_conj_conj_prod_i(complex a,complex b,complex c);
void complex_summ_the_prod(complex a,complex b,complex c);
void complex_summ_the_prod_double(complex a,complex b,double c);
void complex_summ_the_prod_i(complex a,complex b,complex c);
void complex_summ_the_prod_idouble(complex a,complex b,double c);
void complex_summassign(complex a,complex b);
void safe_complex_conj1_prod(complex a,complex b,complex c);
void safe_complex_conj1_prod_minus(complex a,complex b,complex c);
void safe_complex_conj2_prod(complex a,complex b,complex c);
void safe_complex_conj2_prod_minus(complex a,complex b,complex c);
void safe_complex_prod(complex a,complex b,complex c);
void safe_complex_prod_i(complex a,complex b);
void safe_complex_prod_minus(complex a,complex b,complex c);
void safe_complex_prod_minus_i(complex a,complex b);
void unsafe_complex_conj1_prod(complex a,complex b,complex c);
void unsafe_complex_conj1_prod_minus(complex a,complex b,complex c);
void unsafe_complex_conj2_prod(complex a,complex b,complex c);
void unsafe_complex_conj2_prod_minus(complex a,complex b,complex c);
void unsafe_complex_conj_conj_prod(complex a,complex b,complex c);
void unsafe_complex_conj_conj_prod_minus(complex a,complex b,complex c);
void unsafe_complex_prod(complex a,complex b,complex c);
void unsafe_complex_prod_minus(complex a,complex b,complex c);

#endif
