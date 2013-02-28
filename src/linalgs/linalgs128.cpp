#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../new_types/float128.h"
#include "../new_types/new_types_definitions.h"
#include "../routines/mpi.h"


/////////////////////////////////////////////////// quadruple precision ///////////////////////////////////////////

//a=b
void double_vector_from_quadruple_vector(double *a,float_128 *b,int n)
{
#pragma omp for
  for(int i=0;i<n;i++)
    a[i]=double_from_float_128(b[i]);
  
  set_borders_invalid(a);
}

//a=a+b
void quadruple_vector_summassign_double_vector(float_128 *a,double *b,int n)
{
#pragma omp for
  for(int i=0;i<n;i++)
    float_128_summassign_64(a[i],b[i]);

  set_borders_invalid(a);
}

//a=b-c
void quadruple_vector_subt_from_double_vector(float_128 *a,double *b,float_128 *c,int n)
{
#pragma omp for
  for(int i=0;i<n;i++)
    float_128_subt_from_64(a[i],b[i],c[i]);

  set_borders_invalid(a);
}

/////////////////// scalar prodcut in quadruple /////////////////

//(a,b)
void quadruple_vector_glb_scalar_prod(float_128 a,float_128 *b,float_128 *c,int n)
{
#pragma omp single
  {
    float_128 loc_acc={0,0};
    for(int i=0;i<n;i++) float_128_summ_the_prod(loc_acc,b[i],c[i]);
    glb_reduce_float_128(reduce_float_128,loc_acc);
  }
  
  float_128_copy(a,reduce_float_128);
}

//(a,b)
double double_conv_quadruple_vector_glb_scalar_prod(float_128 *a,float_128 *b,int n)
{
  float_128 out;
  quadruple_vector_glb_scalar_prod(out,a,b,n);
  
  return double_from_float_128(out);
}

//////////////// only quadruple accumulation //////////////////

//(a,b)
void quadruple_accumulate_double_vector_glb_scalar_prod(float_128 a,double *b,double *c,int n)
{
#pragma omp single
  {
    float_128 loc_acc={0,0};

    for(int i=0;i<n;i++) float_128_summ_the_64_prod(loc_acc,b[i],c[i]);
    glb_reduce_float_128(reduce_float_128,loc_acc);
  }
  
  float_128_copy(a,reduce_float_128);
}

//(a,b)
double double_conv_quadruple_accumulate_double_vector_glb_scalar_prod(double *a,double *b,int n)
{
  float_128 out;
  quadruple_accumulate_double_vector_glb_scalar_prod(out,a,b,n);
  
  return double_from_float_128(out);
}
