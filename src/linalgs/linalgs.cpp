#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <string.h>
#include <math.h>

#include "../base/communicate.h"
#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../new_types/float128.h"
#include "../new_types/new_types_definitions.h"
#include "../new_types/spin.h"
#include "../new_types/su3.h"
#include "../routines/mpi.h"
#include "../routines/openmp.h"

#ifdef BGP
 #include "../base/bgp_instructions.h"
#endif

//set to zero
void double_vector_init_to_zero(double *a,int n)
{
  if(IS_MASTER_THREAD) memset(a,0,n*sizeof(double));
  
  set_borders_invalid(a);
}

//copy
void double_vector_copy(double *a,double *b,int n)
{
  if(IS_MASTER_THREAD) memcpy(a,b,n*sizeof(double));
  
  set_borders_invalid(a);
}

//summ
THREADABLE_FUNCTION_3ARG(double_vector_summassign, double*,out, double*,in, int,n)
{NISSA_PARALLEL_LOOP(i,n) out[i]+=in[i];set_borders_invalid(out);}}

//prod with double
THREADABLE_FUNCTION_4ARG(double_vector_prod_double, double*,out, double*,in, double,r, int,n)
{NISSA_PARALLEL_LOOP(i,n) out[i]=r*in[i];set_borders_invalid(out);}}

//prod with double of the summ
THREADABLE_FUNCTION_5ARG(double_vector_prod_the_summ_double, double*,out, double,r, double*,in1, double*,in2, int,n)
{NISSA_PARALLEL_LOOP(i,n) out[i]=r*(in1[i]+in2[i]);set_borders_invalid(out);}}

//internal summ
THREADABLE_FUNCTION_4ARG(double_vector_glb_scalar_prod, double*,glb_res, double*,a, double*,b, int,n)
{
#if 0
  //perform thread summ
  double loc_thread_res=0;
  NISSA_PARALLEL_LOOP(i,n)
    loc_thread_res+=a[i]*b[i];
  
  (*glb_res)=glb_reduce_double(loc_thread_res);
#else
  //perform thread summ
  float_128 loc_thread_res={0,0};
  NISSA_PARALLEL_LOOP(i,n)
    float_128_summassign_64(loc_thread_res,a[i]*b[i]);
  
  float_128 temp;
  glb_reduce_float_128(temp,loc_thread_res);
  (*glb_res)=temp[0];
#endif
}}

//put the passed vector to the new norm, returning the reciprocal of normalizating factor
THREADABLE_FUNCTION_5ARG(double_vector_normalize, double*,ratio, double*,out, double*,in, double,norm, int,n)
{
  //compute current norm
  double old_norm;
  double_vector_glb_scalar_prod(&old_norm,in,in,n);
  
  //compute normalizing factor
  double fact=sqrt(norm/old_norm);
  double_vector_prod_double(out,in,fact,n);
  
  if(ratio!=NULL) (*ratio)=1/fact;
}}

//a[]=b[]+c[]*d
THREADABLE_FUNCTION_5ARG(double_vector_summ_double_vector_prod_double, double*,a, double*,b, double*,c, double,d, int,n)
{NISSA_PARALLEL_LOOP(i,n) a[i]=b[i]+c[i]*d;set_borders_invalid(a);}}

//a[]=b[]*c+d[]*e
THREADABLE_FUNCTION_6ARG(double_vector_linear_comb, double*,a, double*,b, double,c, double*,d, double,e, int,n)
{NISSA_PARALLEL_LOOP(i,n) a[i]=b[i]*c+d[i]*e;set_borders_invalid(a);}}

//a[]=b[]-c[]*d
void double_vector_subt_double_vector_prod_double(double *a,double *b,double *c,double d,int n)
{double_vector_summ_double_vector_prod_double(a,b,c,-d,n);}

//////////////////////////////////////////////////////// quadruple precision ///////////////////////////////////////////

//a=b
THREADABLE_FUNCTION_3ARG(double_vector_from_quadruple_vector, double*,a, float_128*,b, int,n)
{NISSA_PARALLEL_LOOP(i,n) a[i]=double_from_float_128(b[i]);set_borders_invalid(a);}}

//a=a+b
THREADABLE_FUNCTION_3ARG(quadruple_vector_summassign_double_vector, float_128*,a, double*,b, int,n)
{NISSA_PARALLEL_LOOP(i,n) float_128_summassign_64(a[i],b[i]);set_borders_invalid(a);}}

//a=b-c
THREADABLE_FUNCTION_4ARG(quadruple_vector_subt_from_double_vector, float_128*,a, double*,b, float_128*,c, int,n)
{NISSA_PARALLEL_LOOP(i,n) float_128_subt_from_64(a[i],b[i],c[i]);set_borders_invalid(a);}}

/////////////////// scalar prodcut in quadruple /////////////////

//(a,b)
void quadruple_vector_glb_scalar_prod(float_128 a,float_128 *b,float_128 *c,int n)
{
  float_128 loc_acc={0,0};
  for(int i=0;i<n;i++) float_128_summ_the_prod(loc_acc,b[i],c[i]);
  glb_reduce_float_128(reduce_float_128,loc_acc);
  
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
  float_128 loc_acc={0,0};
  
  for(int i=0;i<n;i++) float_128_summ_the_64_prod(loc_acc,b[i],c[i]);
  glb_reduce_float_128(reduce_float_128,loc_acc);
  
  float_128_copy(a,reduce_float_128);
}

//(a,b)
double double_conv_quadruple_accumulate_double_vector_glb_scalar_prod(double *a,double *b,int n)
{float_128 out;quadruple_accumulate_double_vector_glb_scalar_prod(out,a,b,n);return double_from_float_128(out);}

//////////////// color put/get from colorspinspin////////////////////////

THREADABLE_FUNCTION_4ARG(get_color_from_colorspinspin, color*,out, colorspinspin*,in, int,id1, int,id2)
{
  NISSA_PARALLEL_LOOP(ivol,loc_vol) get_color_from_colorspinspin(out[ivol],in[ivol],id1,id2);
  set_borders_invalid(out);
}}

THREADABLE_FUNCTION_4ARG(put_color_into_colorspinspin, colorspinspin*,out, color*,in, int,id1, int,id2)
{
  NISSA_PARALLEL_LOOP(ivol,loc_vol) put_color_into_colorspinspin(out[ivol],in[ivol],id1,id2);
  set_borders_invalid(out);
}}

//////////////// color put/get from spincolor//////////////////

THREADABLE_FUNCTION_3ARG(get_color_from_spincolor, color*,out, spincolor*,in, int,id)
{
  NISSA_PARALLEL_LOOP(ivol,loc_vol) get_color_from_spincolor(out[ivol],in[ivol],id);
  set_borders_invalid(out);
}}

THREADABLE_FUNCTION_3ARG(put_color_into_spincolor, spincolor*,out, color*,in, int,id)
{
  NISSA_PARALLEL_LOOP(ivol,loc_vol) put_color_into_spincolor(out[ivol],in[ivol],id);
  set_borders_invalid(out);
}}

//////////////// colorspinspin put/get ////////////////////////

THREADABLE_FUNCTION_3ARG(get_spincolor_from_colorspinspin, spincolor*,out, colorspinspin*,in, int,id)
{
  NISSA_PARALLEL_LOOP(ivol,loc_vol) get_spincolor_from_colorspinspin(out[ivol],in[ivol],id);
  set_borders_invalid(out);
}}

THREADABLE_FUNCTION_3ARG(put_spincolor_into_colorspinspin, colorspinspin*,out, spincolor*,in, int,id)
{
  NISSA_PARALLEL_LOOP(ivol,loc_vol) put_spincolor_into_colorspinspin(out[ivol],in[ivol],id);
  set_borders_invalid(out);
}}

//////////////// su3spinspin put/get ////////////////////////

THREADABLE_FUNCTION_4ARG(get_spincolor_from_su3spinspin, spincolor*,out, su3spinspin*,in, int,id, int,ic)
{
  NISSA_PARALLEL_LOOP(ivol,loc_vol) get_spincolor_from_su3spinspin(out[ivol],in[ivol],id,ic);
  set_borders_invalid(out);
}}

THREADABLE_FUNCTION_4ARG(put_spincolor_into_su3spinspin, su3spinspin*,out, spincolor*,in, int,id, int,ic)
{
  NISSA_PARALLEL_LOOP(ivol,loc_vol) put_spincolor_into_su3spinspin(out[ivol],in[ivol],id,ic);
  set_borders_invalid(out);
}}

////////////////// spincolor algebra/////////////////////

THREADABLE_FUNCTION_3ARG(safe_dirac_prod_spincolor, spincolor*,out, dirac_matr*,m, spincolor*,in)
{
  NISSA_PARALLEL_LOOP(ivol,loc_vol) safe_dirac_prod_spincolor(out[ivol],m,in[ivol]);
  set_borders_invalid(out);
}}
