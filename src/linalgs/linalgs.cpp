#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <math.h>

#include "communicate/communicate.hpp"
#include "base/vectors.hpp"
#include "linalgs/linalgs.hpp"
#include "linalgs/reduce.hpp"
#include "new_types/complex.hpp"
#include "new_types/dirac.hpp"
#include "new_types/float_128.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //set to zero
  void double_vector_init_to_zero(double *a,int64_t n)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      a[i]=0;
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(a);
  }
  
  //copy
  template <class T1,class T2> void internal_vector_copy(T1 *a,T2 *b,int64_t n)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      a[i]=b[i];
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(a);
  }
  
  //double to single and vv
  void double_vector_to_single(float* a,double* b,int64_t n)
  {
    internal_vector_copy(a,b,n);
  }
  
  void single_vector_to_double(double* a,float *b,int64_t n)
  {
    internal_vector_copy(a,b,n);
  }
  
  //double to double
  void double_vector_copy(double *a,double* b,int64_t n)
  {
    internal_vector_copy(a,b,n);
  }
  
  //set to zero
  void single_vector_init_to_zero(float *a,int64_t n)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      a[i]=0;
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(a);
  }
  
  //copy
  void single_vector_copy(float *a,float *b,int64_t n)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      a[i]=b[i];
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(a);
  }
  
  //summ
  void double_vector_summ(double* out,double* in1,double* in2,int64_t n)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      out[i]=in1[i]+in2[i];
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  //subt
  void double_vector_subt(double* out,double* in1,double* in2,int64_t n)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      out[i]=in1[i]-in2[i];
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  //prod with double
  void double_vector_prod_double(double* out,double* in,double r,int64_t n)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      out[i]=r*in[i];
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  void float_128_vector_prod_double(float_128* out,float_128* in,double r,int64_t n)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      float_128_prod_64(out[i],in[i],r);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  //prod with double of the summ
  void double_vector_prod_the_summ_double(double* out,double r,double* in1,double* in2,int64_t n)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      out[i]=r*(in1[i]+in2[i]);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  //scalar product
  void double_vector_glb_scalar_prod(double* glb_res,double* a,double* b,int64_t n)
  {
#ifndef REPRODUCIBLE_RUN
    double *reducing_buffer=get_reducing_buffer<double>(n);
#else
    float_128 *reducing_buffer=get_reducing_buffer<float_128>(n);
#endif
    
    NISSA_PARALLEL_LOOP(i,0,n)
#ifndef REPRODUCIBLE_RUN
      reducing_buffer[i]=a[i]*b[i];
#else
      float_128_from_64(reducing_buffer[i],a[i]*b[i]);
#endif
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
#ifndef REPRODUCIBLE_RUN
    glb_reduce(glb_res,reducing_buffer,n);
#else
    float_128 temp;
    glb_reduce(&temp,reducing_buffer,n);
    *glb_res=temp[0];
#endif
  }
  
  //scalar product
  void complex_vector_glb_scalar_prod(double* glb_res,complex* a,complex* b,int64_t n)
  {
#ifndef REPRODUCIBLE_RUN
    complex *reducing_buffer=get_reducing_buffer<complex>(n);
#else
    complex_128 *reducing_buffer=get_reducing_buffer<complex_128>(n);
#endif
    
    NISSA_PARALLEL_LOOP(i,0,n)
      {
	complex temp;
	unsafe_complex_conj1_prod(temp,a[i],b[i]);
#ifndef REPRODUCIBLE_RUN
	complex_copy(reducing_buffer[i],temp);
#else
	complex_128_from_64(reducing_buffer[i],temp);
#endif
      }
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
#ifndef REPRODUCIBLE_RUN
    glb_reduce((complex*)glb_res,reducing_buffer,n);
#else
    complex_128 temp;
    glb_reduce(&temp,reducing_buffer,n);
    for(int64_t ri=0;ri<2;ri++)
      glb_res[ri]=temp[ri][0];
#endif
  }
  
  //scalar product
  void single_vector_glb_scalar_prod(float* glb_res,float* a,float* b,int64_t n)
  {
    crash("not working");
//     //perform thread summ
//     float loc_thread_res=0;
    
//     NISSA_PARALLEL_LOOP(i,0,n)
//       loc_thread_res+=a[i]*b[i];
//     NISSA_PARALLEL_LOOP_END;
    
//     (*glb_res)=glb_reduce_single(loc_thread_res);
  }
  
  //put a vector of complex equal to its conjugate
  void complex_vector_conj(complex* res,complex* in,int64_t n)
  {
    
    NISSA_PARALLEL_LOOP(i,0,n)
      complex_conj(res[i],in[i]);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(res);
  }
  
  //Ai=Ai+Bi*c
  void complex_vector_summassign_complex_vector_prod_complex(complex* a,complex* b,double* c,int64_t n)
  {
    
    NISSA_PARALLEL_LOOP(i,0,n)
      complex_summ_the_prod(a[i],b[i],c);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(a);
  }
  
  //put the passed vector to the new norm, returning the reciprocal of normalizating factor
  void double_vector_normalize(double* ratio,double* out,double* in,double norm,int64_t n)
  {
    //compute current norm
    double old_norm;
    double_vector_glb_scalar_prod(&old_norm,in,in,n);
    
    //compute normalizing factor
    double fact=sqrt(norm/old_norm);
    double_vector_prod_double(out,in,fact,n);
    
    if(ratio!=NULL) (*ratio)=1/fact;
    
    set_borders_invalid(out);
  }
  
  //a[]=b[]+c[]*d
  void double_vector_summ_double_vector_prod_double(double* a,double* b,double* c,double d,int64_t n,int OPT)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      a[i]=b[i]+c[i]*d;
    NISSA_PARALLEL_LOOP_END;
    
    if(!(OPT&DO_NOT_SET_FLAGS)) set_borders_invalid(a);
  }
  
  //single version
  void single_vector_summ_single_vector_prod_single(float* a,float* b,float* c,float d,int64_t n,int OPT)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      a[i]=b[i]+c[i]*d;
    NISSA_PARALLEL_LOOP_END;
    
    if(!(OPT&DO_NOT_SET_FLAGS)) set_borders_invalid(a);
  }
  
  //a[]=b[]*c+d[]*e
  void double_vector_linear_comb(double* a,double* b,double c,double* d,double e,int64_t n,int OPT)
  {
    
    NISSA_PARALLEL_LOOP(i,0,n)
      a[i]=b[i]*c+d[i]*e;
    NISSA_PARALLEL_LOOP_END;
    
    if(!(OPT&DO_NOT_SET_FLAGS))
      set_borders_invalid(a);
  }
  
  //a[]=b[]*c+d[]*e
  void single_vector_linear_comb(float* a,float* b,float c,float* d,float e,int64_t n,int OPT)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      a[i]=b[i]*c+d[i]*e;
    NISSA_PARALLEL_LOOP_END;
    
    if(!(OPT&DO_NOT_SET_FLAGS)) set_borders_invalid(a);
  }
  
  ////////////////////////////////////////////////////// quadruple precision ///////////////////////////////////////////
  
  //a=b
  void double_vector_from_quadruple_vector(double* a,float_128* b,int64_t n)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      a[i]=double_from_float_128(b[i]);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(a);
  }
  void quadruple_vector_from_double_vector(float_128* a,double* b,int64_t n)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      float_128_from_64(a[i],b[i]);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(a);
  }
  
  //a=a+b
  void quadruple_vector_summassign_double_vector(float_128* a,double* b,int64_t n)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      float_128_summassign_64(a[i],b[i]);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(a);
  }
  
  //a=b-c
  void quadruple_vector_subt_from_double_vector(float_128* a,double* b,float_128* c,int64_t n)
  {
    NISSA_PARALLEL_LOOP(i,0,n)
      float_128_subt_from_64(a[i],b[i],c[i]);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(a);
  }
  
  /////////////////// scalar prodcut in quadruple /////////////////
  
  //(a,b)
  void quadruple_vector_glb_scalar_prod(float_128* glb_res,float_128* a,float_128* b,int64_t n)
  {
    crash("");
    glb_res[0][0]=0; //to remove warning
    // //perform thread summ
    // float_128 loc_thread_res={0,0};
    // NISSA_PARALLEL_LOOP(i,0,n)
    //   float_128_summ_the_prod(loc_thread_res,a[i],b[i]);
    // NISSA_PARALLEL_LOOP_END;
    
    // glb_reduce_float_128(*glb_res,loc_thread_res);
  }
  
  //(a,b)
  void double_conv_quadruple_vector_glb_scalar_prod(double *out,float_128 *a,float_128 *b,int64_t n)
  {
    float_128 out_128;
    quadruple_vector_glb_scalar_prod(&out_128,a,b,n);
    *out=out_128[0];
  }
  
  //////////////// only quadruple accumulation //////////////////
  
  //(a,b)
  void quadruple_accumulate_double_vector_glb_scalar_prod(float_128* a,double* b,double* c,int64_t n)
  {
    crash("");
    a[0][0]=a[0][1]=0; //to remove warning
    // //perform thread summ
    // float_128 loc_thread_res={0,0};
    // NISSA_PARALLEL_LOOP(i,0,n)
    //   float_128_summassign_64(loc_thread_res,b[i]*c[i]);
    // NISSA_PARALLEL_LOOP_END;
    
    // glb_reduce_float_128(*a,loc_thread_res);
  }
  
  //(a,b)
  double double_conv_quadruple_accumulate_double_vector_glb_scalar_prod(double *a,double *b,int64_t n)
  {
    float_128 out;
    quadruple_accumulate_double_vector_glb_scalar_prod(&out,a,b,n);
    return double_from_float_128(out);
  }
  
  ////////////////////// color put/get from su3 ///////////////////////////
  
  void get_color_from_su3(eo_ptr<color> out,eo_ptr<su3> in,int ic_source)
  {
    for(int eo=0;eo<2;eo++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,locVolh)
	  get_color_from_su3(out[eo][ivol],in[eo][ivol],ic_source);
	NISSA_PARALLEL_LOOP_END;
	set_borders_invalid(out[eo]);
      }
  }
  
    void put_color_into_su3(eo_ptr<su3> out,eo_ptr<color> in,int ic_source)
  {
    for(int eo=0;eo<2;eo++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,locVolh)
	  put_color_into_su3(out[eo][ivol],in[eo][ivol],ic_source);
	NISSA_PARALLEL_LOOP_END;
	set_borders_invalid(out[eo]);
      }
  }
  
  //////////////// color put/get from colorspinspin////////////////////////
  
  void get_color_from_colorspinspin(color* out,colorspinspin* in,int id1,int id2)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      get_color_from_colorspinspin(out[ivol],in[ivol],id1,id2);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  void put_color_into_colorspinspin(colorspinspin* out,color* in,int id1,int id2)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      put_color_into_colorspinspin(out[ivol],in[ivol],id1,id2);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  //////////////// color put/get from spincolor//////////////////
  
  void get_color_from_spincolor(color* out,spincolor* in,int id)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      get_color_from_spincolor(out[ivol],in[ivol],id);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  void put_color_into_spincolor(spincolor* out,color* in,int id)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      put_color_into_spincolor(out[ivol],in[ivol],id);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  //////////////// colorspinspin put/get ////////////////////////
  
  void get_spincolor_from_colorspinspin(spincolor* out,colorspinspin* in,int id)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      get_spincolor_from_colorspinspin(out[ivol],in[ivol],id);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  void put_spincolor_into_colorspinspin(colorspinspin* out,spincolor* in,int id)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      put_spincolor_into_colorspinspin(out[ivol],in[ivol],id);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  //////////////// su3spinspin put/get ////////////////////////
  
  void get_spincolor_from_su3spinspin(spincolor* out,su3spinspin* in,int id,int ic)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      get_spincolor_from_su3spinspin(out[ivol],in[ivol],id,ic);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  void put_spincolor_into_su3spinspin(su3spinspin* out,spincolor* in,int id,int ic)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      put_spincolor_into_su3spinspin(out[ivol],in[ivol],id,ic);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  ////////////////// spincolor algebra/////////////////////
  
  void safe_dirac_prod_spincolor(spincolor* out,const dirac_matr& m,spincolor* in)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      safe_dirac_prod_spincolor(out[ivol],m,in[ivol]);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  void safe_dirac_prod_colorspinspin(colorspinspin* out,const dirac_matr& m,colorspinspin* in)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      safe_dirac_prod_colorspinspin(out[ivol],m,in[ivol]);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  ///////////////////// rotations ////////////////////////
  
  void rotate_vol_colorspinspin_to_physical_basis(colorspinspin* s,int rsi,int rso)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int ic=0;ic<3;ic++)
	rotate_spinspin_to_physical_basis(s[ivol][ic],rsi,rso);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(s);
  }
  
  void rotate_vol_su3spinspin_to_physical_basis(su3spinspin* s,int rsi,int rso)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int ic1=0;ic1<3;ic1++)
	for(int ic2=0;ic2<3;ic2++)
	  rotate_spinspin_to_physical_basis(s[ivol][ic1][ic2],rsi,rso);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(s);
  }
  
  //ildg to nissa and vice-versa
  void quad_su3_nissa_to_ildg_reord_in_place(quad_su3* in)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      quad_su3_nissa_to_ildg_reord(in[ivol],in[ivol]);
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
  }
  void quad_su3_ildg_to_nissa_reord_in_place(quad_su3* in)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      quad_su3_ildg_to_nissa_reord(in[ivol],in[ivol]);
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
  }
  
  void parallel_memcpy(void* out,void* in,int64_t n)
  {
#if THREADS_TYPE == OPENMP_THREADS
    
#pragma omp parallel
    {
      NISSA_CHUNK_WORKLOAD(start,chunk_load,end,0,n,THREAD_ID,nthreads);
      memcpy((char*)out+start,(char*)in+start,chunk_load);
      (void)end;//to avoid warning
    }
    THREAD_BARRIER();
#else
    memcpy(out,in,n);
#endif
  }
}
