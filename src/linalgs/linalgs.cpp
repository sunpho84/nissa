#include <string.h>

#include "../new_types/new_types_definitions.h"
#include "../new_types/su3.h"
#include "../new_types/spin.h"
#include "../new_types/float128.h"
#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../base/routines.h"
#include "../base/communicate.h"

#ifdef BGP
 #include "../base/bgp_instructions.h"
#endif

//set to zero
void double_vector_init_to_zero(double *a,int n)
{
#pragma omp single
  memset(a,0,n*sizeof(double));
  
  set_borders_invalid(a);
}

//copy
void double_vector_copy(double *a,double *b,int n)
{
#pragma omp single
  memcpy(a,b,n*sizeof(double));
  
  set_borders_invalid(a);
}

//summ
void double_vector_summassign(double *out,double *in,int n)
{
#pragma omp for
  for(int i=0;i<n;i++)
    out[i]+=in[i];
  
  set_borders_invalid(out);
}

//prod with double
void double_vector_prod_double(double *out,double *in,double r,int n)
{
#pragma omp for
  for(int i=0;i<n;i++)
    out[i]=r*in[i];
  
  set_borders_invalid(out);
}

//prod with double of the summ
void double_vector_prod_the_summ_double(double *out,double r,double *in1,double *in2,int n)
{
#pragma omp for
  for(int i=0;i<n;i++)
    out[i]=r*(in1[i]+in2[i]);
  
  set_borders_invalid(out);
}

//scalar product between a and b
double double_vector_loc_scalar_prod(double *a,double *b,int n)
{
#pragma omp single
  reduce_double=0;
#pragma omp for reduction(+:reduce_double)
  for(int i=0;i<n;i++)
    reduce_double+=a[i]*b[i];
  
  return reduce_double;
}

//a[]=b[]+c[]*d
void double_vector_summ_double_vector_prod_double(double *a,double *b,double *c,double d,int n)
{
#pragma omp for
  for(int i=0;i<n;i++)
    a[i]=b[i]+c[i]*d;
  
  set_borders_invalid(a);
}

//a[]=b[]*c+d[]*e
void double_vector_linear_comb(double *a,double *b,double c,double *d,double e,int n)
{
#pragma omp for
  for(int i=0;i<n;i++)
    a[i]=b[i]*c+d[i]*e;
  
  set_borders_invalid(a);
}

//a[]=b[]-c[]*d
void double_vector_subt_double_vector_prod_double(double *a,double *b,double *c,double d,int n)
{double_vector_summ_double_vector_prod_double(a,b,c,-d,n);}

//global reduction wrapper
double double_vector_glb_scalar_prod(double *a,double *b,int n)
{return glb_reduce_double(double_vector_loc_scalar_prod(a,b,n));}

//////////////////////////////////////////////////////// quadruple precision ///////////////////////////////////////////

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

//////////////// color put/get from colorspinspin////////////////////////

void get_color_from_colorspinspin(color *out,colorspinspin *in,int id1,int id2)
{
  #pragma omp for
  nissa_loc_vol_loop(ivol)
    get_color_from_colorspinspin(out[ivol],in[ivol],id1,id2);
  
  set_borders_invalid(out);
}

void put_color_into_colorspinspin(colorspinspin *out,color *in,int id1,int id2)
{
#pragma omp for
  nissa_loc_vol_loop(ivol)
    put_color_into_colorspinspin(out[ivol],in[ivol],id1,id2);
  
  set_borders_invalid(out);
}

//////////////// color put/get from spincolor////////////////////////

void get_color_from_spincolor(color *out,spincolor *in,int id)
{
#pragma omp for
  nissa_loc_vol_loop(ivol)
    get_color_from_spincolor(out[ivol],in[ivol],id);
  
  set_borders_invalid(out);
}

void put_color_into_spincolor(spincolor *out,color *in,int id)
{
#pragma omp for
  nissa_loc_vol_loop(ivol)
    put_color_into_spincolor(out[ivol],in[ivol],id);
  
  set_borders_invalid(out);
}

//////////////// colorspinspin put/get ////////////////////////

void get_spincolor_from_colorspinspin(spincolor *out,colorspinspin *in,int id)
{
#pragma omp for
  nissa_loc_vol_loop(ivol)
    get_spincolor_from_colorspinspin(out[ivol],in[ivol],id);
  
  set_borders_invalid(out);
}

void put_spincolor_into_colorspinspin(colorspinspin *out,spincolor *in,int id)
{
#pragma omp for
  nissa_loc_vol_loop(ivol)
    put_spincolor_into_colorspinspin(out[ivol],in[ivol],id);
  
  set_borders_invalid(out);
}

//////////////// su3spinspin put/get ////////////////////////

void get_spincolor_from_su3spinspin(spincolor *out,su3spinspin *in,int id,int ic)
{
#pragma omp for
  nissa_loc_vol_loop(ivol)
    get_spincolor_from_su3spinspin(out[ivol],in[ivol],id,ic);
  
  set_borders_invalid(out);
}

void put_spincolor_into_su3spinspin(su3spinspin *out,spincolor *in,int id,int ic)
{
#pragma omp for
  nissa_loc_vol_loop(ivol)
    put_spincolor_into_su3spinspin(out[ivol],in[ivol],id,ic);
  
  set_borders_invalid(out);
}

////////////////// spincolor algebra/////////////////////

void safe_dirac_prod_spincolor(spincolor *out,dirac_matr &m,spincolor *in)
{
#pragma omp for
  nissa_loc_vol_loop(ivol)
    safe_dirac_prod_spincolor(out[ivol],&m,in[ivol]);
  
  set_borders_invalid(out);
}
