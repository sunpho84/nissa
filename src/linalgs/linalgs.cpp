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
  memset(a,0,n*sizeof(double));
  set_borders_invalid(a);
}

//copy
void double_vector_copy(double *a,double *b,int n)
{
  memcpy(a,b,n*sizeof(double));
  set_borders_invalid(a);
}

//scalar product between a and b
double double_vector_loc_scalar_prod(double *a,double *b,int n)
{
  double out=0;
  
  int i=0;
  
#ifdef BGP
  bgp_complex N0,N1,N2;
  bgp_color_put_to_zero(N0,N1,N2);
  while(i+24<n)
    {
      bgp_complex A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,AA,AB;
      bgp_complex B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,BA,BB;
      bgp_spincolor_load(A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,AA,AB,(*((spincolor*)(a+i))));
      bgp_spincolor_load(B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,BA,BB,(*((spincolor*)(b+i))));
      bgp_summassign_color_realscalarprod_spincolor(N0,N1,N2, A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,AA,AB, B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,BA,BB);
      
      i+=24;
    }
  bgp_square_norm_color(N0,N1,N2);
  complex cout;
  bgp_complex_save(cout,N0);
  out=cout[0];
#endif
  
  //reduce the remaining data (or all if not on bgp)
  while(i<n)
    {
      out+=a[i]*b[i];
      i++;
    }
  
  return out;
}

//a[]=b[]+c[]*d
void double_vector_summ_double_vector_prod_double(double *a,double *b,double *c,double d,int n)
{
  for(int i=0;i<n;i++) a[i]=b[i]+c[i]*d;
  set_borders_invalid(a);
}

//a[]=b[]*c+d[]*e
void double_vector_linear_comb(double *a,double *b,double c,double *d,double e,int n)
{
  for(int i=0;i<n;i++) a[i]=b[i]*c+d[i]*e;
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
  for(int i=0;i<n;i++) a[i]=double_from_float_128(b[i]);

  set_borders_invalid(a);
}

//a=a+b
void quadruple_vector_summassign_double_vector(float_128 *a,double *b,int n)
{
  for(int i=0;i<n;i++) float_128_summassign_64(a[i],b[i]);

  set_borders_invalid(a);
}

//a=b-c
void quadruple_vector_subt_from_double_vector(float_128 *a,double *b,float_128 *c,int n)
{
  for(int i=0;i<n;i++) float_128_subt_from_64(a[i],b[i],c[i]);

  set_borders_invalid(a);
}

/////////////////// scalar prodcut in quadruple /////////////////

//(a,b)
void quadruple_vector_glb_scalar_prod(float_128 a,float_128 *b,float_128 *c,int n)
{
  float_128 loc_a={0,0};
  for(int i=0;i<n;i++) float_128_summ_the_prod(loc_a,b[i],c[i]);
  MPI_Allreduce(loc_a,a,1,MPI_FLOAT_128,MPI_FLOAT_128_SUM,MPI_COMM_WORLD);
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
  float_128 loc_a={0,0};
  for(int i=0;i<n;i++) float_128_summ_the_64_prod(loc_a,b[i],c[i]);
  MPI_Allreduce(loc_a,a,1,MPI_FLOAT_128,MPI_FLOAT_128_SUM,MPI_COMM_WORLD);
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
  nissa_loc_vol_loop(ivol)
    get_color_from_colorspinspin(out[ivol],in[ivol],id1,id2);
  set_borders_invalid(out);
}

void put_color_into_colorspinspin(colorspinspin *out,color *in,int id1,int id2)
{
  nissa_loc_vol_loop(ivol)
    put_color_into_colorspinspin(out[ivol],in[ivol],id1,id2);
  set_borders_invalid(out);
}

//////////////// color put/get from spincolor////////////////////////

void get_color_from_spincolor(color *out,spincolor *in,int id)
{
  nissa_loc_vol_loop(ivol)
    get_color_from_spincolor(out[ivol],in[ivol],id);
  set_borders_invalid(out);
}

void put_color_into_spincolor(spincolor *out,color *in,int id)
{
  nissa_loc_vol_loop(ivol)
    put_color_into_spincolor(out[ivol],in[ivol],id);
  set_borders_invalid(out);
}

//////////////// colorspinspin put/get ////////////////////////

void get_spincolor_from_colorspinspin(spincolor *out,colorspinspin *in,int id)
{
  nissa_loc_vol_loop(ivol)
    get_spincolor_from_colorspinspin(out[ivol],in[ivol],id);
  set_borders_invalid(out);
}

void put_spincolor_into_colorspinspin(colorspinspin *out,spincolor *in,int id)
{
  nissa_loc_vol_loop(ivol)
    put_spincolor_into_colorspinspin(out[ivol],in[ivol],id);
  set_borders_invalid(out);
}

//////////////// su3spinspin put/get ////////////////////////

void get_spincolor_from_su3spinspin(spincolor *out,su3spinspin *in,int id,int ic)
{
  nissa_loc_vol_loop(ivol)
    get_spincolor_from_su3spinspin(out[ivol],in[ivol],id,ic);
  set_borders_invalid(out);
}

void put_spincolor_into_su3spinspin(su3spinspin *out,spincolor *in,int id,int ic)
{
  nissa_loc_vol_loop(ivol)
    put_spincolor_into_su3spinspin(out[ivol],in[ivol],id,ic);
  set_borders_invalid(out);
}

////////////////// spincolor algebra/////////////////////

void safe_dirac_prod_spincolor(spincolor *out,dirac_matr &m,spincolor *in)
{
  nissa_loc_vol_loop(ivol)
    safe_dirac_prod_spincolor(out[ivol],&m,in[ivol]);
  set_borders_invalid(out);
}
