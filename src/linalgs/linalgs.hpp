#ifndef _LINALGS_HPP
#define _LINALGS_HPP

#include "new_types/dirac.hpp"
#include "new_types/float_128.hpp"
#include "new_types/su3.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  void complex_vector_glb_collapse(double *res,complex *a,int n);
  void complex_vector_conj(complex *res,complex *in,int n);
  inline void complex_vector_self_conj(complex *v,int n)
  {complex_vector_conj(v,v,n);}
  void complex_vector_glb_scalar_prod(double *glb_res,complex *a,complex *b,int n);
  void complex_vector_summassign_complex_vector_prod_complex(complex *a,complex *b,complex c,int n);
  inline void complex_vector_subtassign_complex_vector_prod_complex(complex *a,complex *b,complex c,int n)
  {
    complex d={-c[RE],-c[IM]};
    complex_vector_summassign_complex_vector_prod_complex(a,b,d,n);
  }
  
  void single_vector_init_to_zero(float *a,int n);
  void single_vector_copy(float *a,float *b,int n);
  void double_conv_quadruple_vector_glb_scalar_prod(double *out,float_128 *a,float_128 *b,int n);
  void double_vector_glb_scalar_prod(double *res,double *a,double *b,int n);
  
  template <class T> double double_vector_glb_norm2(T *v,int n_per_class)
  {double res;double_vector_glb_scalar_prod(&res,(double*)v,(double*)v,n_per_class*sizeof(T)/sizeof(double));return res;}
  
  //Norm2 of a scalar
  template <typename T>
  T norm2(T &s)
  {
    return s*s;
  }
  
  //Norm2 of an array
  template <typename TIn,
	    int N,
	    typename TOut=std::remove_all_extents_t<TIn>>
  TOut norm2(TIn (&v)[N])
  {
    TOut out=0;
    
    for(int i=0;i<N;i++)
      out+=norm2(v[i]);
    
    return out;
  }
  
  //Takes the trace of the square of the vector, on all internal T indices
  template <typename TOut,
	    typename TIn,
	    typename=std::enable_if_t<std::is_same<TOut,typename std::remove_all_extents<TIn>::type>::value>>
  void vector_loc_norm2(TOut *loc_norm2,TIn *v,int n_per_class)
  {
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(i,0,n_per_class)
      loc_norm2[i]+=norm2(v[i]);
    
    set_borders_invalid(loc_norm2);
  }
  
  void single_vector_glb_scalar_prod(float *res,float *a,float *b,int n);
  void double_vector_glb_collapse(double *res,double *a,int n);
  void double_vector_copy(double *a,double *b,int n);
  void double_vector_to_single(float *a,double *b,int n);
  void single_vector_to_double(double *a,float *b,int n);
  void double_vector_from_quadruple_vector(double *a,float_128 *b,int n);
  void quadruple_vector_from_double_vector(float_128 *a,double *b,int n);
  void double_vector_init_to_zero(double *a,int n);
  void double_vector_linear_comb(double *a,double *b,double c,double *d,double e,int n,int OPT=0);
  void single_vector_linear_comb(float *a,float *b,float c,float *d,float e,int n,int OPT=0);
  void double_vector_prod_double(double *out,double *in,double r,int n);
  inline void double_vector_prodassign_double(double *v,double r,int n)
  {double_vector_prod_double(v,v,r,n);}
  void double_vector_normalize(double *ratio,double *out,double *in,double norm,int n);
  void double_vector_prod_the_summ_double(double *out,double r,double *in1,double *in2,int n);
  void double_vector_subt(double *out,double *in1,double *i2,int n);
  void double_vector_summ(double *out,double *in1,double *i2,int n);
  inline void double_vector_summassign(double *out,double *in,int n){double_vector_summ(out,out,in,n);}
  inline void double_vector_subtassign(double *out,double *in,int n){double_vector_subt(out,out,in,n);}
  void double_vector_summ_double_vector_prod_double(double *a,double *b,double *c,double d,int n,int OPT=0);
  inline void double_vector_summassign_double_vector_prod_double(double *a,double *b,double c,int n,int OPT=0)
  {double_vector_summ_double_vector_prod_double(a,a,b,c,n,OPT);}
  void float_128_vector_prod_double(float_128 *out,float_128 *in,double r,int n);
  void single_vector_summ_single_vector_prod_single(float *a,float *b,float *c,float d,int n,int OPT=0);
  void get_color_from_colorspinspin(color *out,colorspinspin *in,int id1,int id2);
  void get_color_from_spincolor(color *out,spincolor *in,int id);
  void get_color_from_su3(color **out,su3 **in,int ic);
  void get_spincolor_from_colorspinspin(spincolor *out,colorspinspin *in,int id);
  void get_spincolor_from_su3spinspin(spincolor *out,su3spinspin *in,int id,int ic);
  void parallel_memcpy(void *out,void *in,int n);
  void put_color_into_su3(su3 **out,color **in,int ic);
  void put_color_into_colorspinspin(colorspinspin *out,color *in,int id1,int id2);
  void put_color_into_spincolor(spincolor *out,color *in,int id);
  void put_spincolor_into_colorspinspin(colorspinspin *out,spincolor *in,int id);
  void put_spincolor_into_su3spinspin(su3spinspin *out,spincolor *in,int id,int ic);
  void quadruple_accumulate_double_vector_glb_scalar_prod(float_128 a,double *b,double *c,int n);
  void quadruple_vector_glb_scalar_prod(float_128 a,float_128 *b,float_128 *c,int n);
  void quadruple_vector_subt_from_double_vector(float_128 *a,double *b,float_128 *c,int n);
  void quadruple_vector_summassign_double_vector(float_128 *a,double *b,int n);
  void safe_dirac_prod_spincolor(spincolor *out,dirac_matr *m,spincolor *in);
  void safe_dirac_prod_colorspinspin(colorspinspin *out,dirac_matr *m,colorspinspin *in);
  void rotate_vol_colorspinspin_to_physical_basis(colorspinspin *s,int rsi,int rso);
  void rotate_vol_su3spinspin_to_physical_basis(su3spinspin *s,int rsi,int rso);
  void quad_su3_nissa_to_ildg_reord_in_place(quad_su3 *conf);
  void quad_su3_ildg_to_nissa_reord_in_place(quad_su3 *conf);
}

#endif
