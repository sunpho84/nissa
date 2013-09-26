#ifndef _LINALGS_H
#define _LINALGS_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void double_conv_quadruple_vector_glb_scalar_prod(double *out,float_128 *a,float_128 *b,int n);
  void double_vector_glb_scalar_prod(double *res,double *a,double *b,int n);
  void complex_vector_glb_collapse(double *res,complex *a,int n);
  void double_vector_glb_collapse(double *res,double *a,int n);
  void double_vector_copy(double *a,double *b,int n);
  void double_vector_from_quadruple_vector(double *a,float_128 *b,int n);
  void quadruple_vector_from_double_vector(float_128 *a,double *b,int n);
  void double_vector_init_to_zero(double *a,int n);
  void double_vector_linear_comb(double *a,double *b,double c,double *d,double e,int n,int OPT=0);
  void double_vector_prod_double(double *out,double *in,double r,int n);
  void double_vector_normalize(double *rat,double *out,double *in,double fact,int n);
  void double_vector_prod_the_summ_double(double *out,double r,double *in1,double *in2,int n);
  void double_vector_summassign(double *out,double *in,int n);
  void double_vector_subt(double *out,double *in1,double *i2,int n);
  void double_vector_summ_double_vector_prod_double(double *a,double *b,double *c,double d,int n,int OPT=0);
  void get_color_from_colorspinspin(color *out,colorspinspin *in,int id1,int id2);
  void get_color_from_spincolor(color *out,spincolor *in,int id);
  void get_spincolor_from_colorspinspin(spincolor *out,colorspinspin *in,int id);
  void get_spincolor_from_su3spinspin(spincolor *out,su3spinspin *in,int id,int ic);
  void parallel_memcpy(void *out,void *in,int n);
  void put_color_into_colorspinspin(colorspinspin *out,color *in,int id1,int id2);
  void put_color_into_spincolor(spincolor *out,color *in,int id);
  void put_spincolor_into_colorspinspin(colorspinspin *out,spincolor *in,int id);
  void put_spincolor_into_su3spinspin(su3spinspin *out,spincolor *in,int id,int ic);
  void quadruple_accumulate_double_vector_glb_scalar_prod(float_128 a,double *b,double *c,int n);
  void quadruple_vector_glb_scalar_prod(float_128 a,float_128 *b,float_128 *c,int n);
  void quadruple_vector_subt_from_double_vector(float_128 *a,double *b,float_128 *c,int n);
  void quadruple_vector_summassign_double_vector(float_128 *a,double *b,int n);
  void safe_dirac_prod_spincolor(spincolor *out,dirac_matr *m,spincolor *in);
  void rotate_vol_colorspinspin_to_physical_basis(colorspinspin *s,int rsi,int rso);
  void rotate_vol_su3spinspin_to_physical_basis(su3spinspin *s,int rsi,int rso);
}

#endif
