#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <math.h>

#include "communicate/communicate.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "new_types/complex.hpp"
#include "new_types/float128.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif
#ifdef BGQ
 #include "bgq/intrinsic.hpp"
#endif

namespace nissa
{
  //set to zero
  void double_vector_init_to_zero(double *a,int n)
  {
    GET_THREAD_ID();
    if(IS_MASTER_THREAD) memset(a,0,n*sizeof(double));
    
    set_borders_invalid(a);
  }
  
  //copy
  void double_vector_copy(double *a,double *b,int n)
  {
    GET_THREAD_ID();
    if(IS_MASTER_THREAD) memcpy(a,b,n*sizeof(double));
    
    set_borders_invalid(a);
  }
  
  //summ
  THREADABLE_FUNCTION_3ARG(double_vector_summassign, double*,out, double*,in, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) out[i]+=in[i];set_borders_invalid(out);}}

  //subt
  THREADABLE_FUNCTION_4ARG(double_vector_subt, double*,out, double*,in1, double*,in2, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) out[i]=in1[i]-in2[i];set_borders_invalid(out);}}
  
  //prod with double
  THREADABLE_FUNCTION_4ARG(double_vector_prod_double, double*,out, double*,in, double,r, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) out[i]=r*in[i];set_borders_invalid(out);}}

  //prod with double of the summ
  THREADABLE_FUNCTION_5ARG(double_vector_prod_the_summ_double, double*,out, double,r, double*,in1, double*,in2, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) out[i]=r*(in1[i]+in2[i]);set_borders_invalid(out);}}
  
  //scalar product
  THREADABLE_FUNCTION_4ARG(double_vector_glb_scalar_prod, double*,glb_res, double*,a, double*,b, int,n)
  {
#ifndef REPRODUCIBLE_RUN
    //perform thread summ
    double loc_thread_res=0;
    GET_THREAD_ID();
    
#ifndef BGQ
    NISSA_PARALLEL_LOOP(i,0,n)
      loc_thread_res+=a[i]*b[i];
#else
    int max_n=n/8;
    DECLARE_REG_BI_HALFSPIN(reg_loc_thread_res);
    REG_SPLAT_BI_HALFSPIN(reg_loc_thread_res,0);
    NISSA_PARALLEL_LOOP(i,0,max_n)
      {
	DECLARE_REG_BI_HALFSPIN(reg_a);
	DECLARE_REG_BI_HALFSPIN(reg_b);
	
	double *a_ptr=(double*)(((bi_halfspin*)a)[i]);
	double *b_ptr=(double*)(((bi_halfspin*)b)[i]);
	
	BI_HALFSPIN_PREFETCH_NEXT_NEXT(a);
	BI_HALFSPIN_PREFETCH_NEXT_NEXT(b);
	REG_LOAD_BI_HALFSPIN(reg_a,a_ptr);
	REG_LOAD_BI_HALFSPIN(reg_b,b_ptr);
	REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(reg_loc_thread_res,s0),NAME2(reg_loc_thread_res,s0),
					     NAME2(reg_a,s0),NAME2(reg_b,s0));
	REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(reg_loc_thread_res,s1),NAME2(reg_loc_thread_res,s1),
					     NAME2(reg_a,s1),NAME2(reg_b,s1));
      }
    REG_BI_COMPLEX_SUMM(reg_loc_thread_res_s0,reg_loc_thread_res_s0,reg_loc_thread_res_s1);
    
    //store
    bi_complex temp;
    STORE_REG_BI_COMPLEX(temp,reg_loc_thread_res_s0);
    loc_thread_res=temp[0][0]+temp[0][1]+temp[1][0]+temp[1][1];
    
    //last part
    if(IS_MASTER_THREAD) for(int i=max_n*8;i<n;i++) loc_thread_res+=a[i]*b[i];
#endif
    
    (*glb_res)=glb_reduce_double(loc_thread_res);
    
#else //reproducible run
    
    //perform thread summ
    float_128 loc_thread_res={0,0};
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(i,0,n)
      float_128_summassign_64(loc_thread_res,a[i]*b[i]);
    
    float_128 temp;
    glb_reduce_float_128(temp,loc_thread_res);
    (*glb_res)=temp[0];
#endif
  }}

  //summ all points
  THREADABLE_FUNCTION_3ARG(double_vector_glb_collapse, double*,glb_res, double*,a, int,n)
  {
#ifndef REPRODUCIBLE_RUN
    //perform thread summ
    double loc_thread_res=0;
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(i,0,n)
      loc_thread_res+=a[i];
    
    (*glb_res)=glb_reduce_double(loc_thread_res);
#else
    //perform thread summ
    float_128 loc_thread_res={0,0};
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(i,0,n)
      float_128_summassign_64(loc_thread_res,a[i]);
    
    float_128 temp;
    glb_reduce_float_128(temp,loc_thread_res);
    (*glb_res)=temp[0];
#endif
  }}

  //complex version
  THREADABLE_FUNCTION_3ARG(complex_vector_glb_collapse, double*,glb_res, complex*,a, int,n)
  {
#ifndef REPRODUCIBLE_RUN
    //perform thread summ
    complex loc_thread_res={0,0};
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(i,0,n)
      complex_summassign(loc_thread_res,a[i]);
    
    for(int ri=0;ri<2;ri++)
      glb_res[ri]=glb_reduce_double(loc_thread_res[ri]);
#else
    //perform thread summ
    complex_128 loc_thread_res={{0,0},{0,0}};
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(i,0,n)
      complex_128_summassign_64(loc_thread_res,a[i]);
    
    //drop back to complex after reducing all threads and ranks
    for(int ri=0;ri<2;ri++)
      {
	float_128 temp;
	glb_reduce_float_128(temp,loc_thread_res[ri]);
	glb_res[ri]=temp[0];
      }
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
    
    set_borders_invalid(out);
  }}

  //a[]=b[]+c[]*d
  THREADABLE_FUNCTION_6ARG(double_vector_summ_double_vector_prod_double, double*,a, double*,b, double*,c, double,d, int,n, int,OPT)
  {
    GET_THREAD_ID();
#ifndef BGQ
    NISSA_PARALLEL_LOOP(i,0,n) a[i]=b[i]+c[i]*d;
#else
    int max_n=n/8;
    DECLARE_REG_BI_COMPLEX(reg_d);
    REG_SPLAT_BI_COMPLEX(reg_d,d);
    NISSA_PARALLEL_LOOP(i,0,max_n)
      {
	DECLARE_REG_BI_HALFSPIN(reg_in1);
	DECLARE_REG_BI_HALFSPIN(reg_in2);
	DECLARE_REG_BI_HALFSPIN(reg_out);
	
	double *in1=(double*)(((bi_halfspin*)b)[i]);
	double *in2=(double*)(((bi_halfspin*)c)[i]);
	
	BI_HALFSPIN_PREFETCH_NEXT_NEXT(in1);
	BI_HALFSPIN_PREFETCH_NEXT_NEXT(in2);
	
	REG_LOAD_BI_HALFSPIN(reg_in1,in1);
	REG_LOAD_BI_HALFSPIN(reg_in2,in2);
	REG_BI_HALFSPIN_SUMM_THE_PROD_4DOUBLE(reg_out,reg_in1,reg_in2,reg_d);
	STORE_REG_BI_HALFSPIN(((bi_halfspin*)a)[i],reg_out);
      }
    //last part
    if(IS_MASTER_THREAD) for(int i=max_n*8;i<n;i++) a[i]=b[i]+c[i]*d;
#endif
    if(!(OPT&DO_NOT_SET_FLAGS)) set_borders_invalid(a);
  }}

  //a[]=b[]*c+d[]*e
  THREADABLE_FUNCTION_7ARG(double_vector_linear_comb, double*,a, double*,b, double,c, double*,d, double,e, int,n, int,OPT)
  {
    GET_THREAD_ID();
#ifndef BGQ
    NISSA_PARALLEL_LOOP(i,0,n) a[i]=b[i]*c+d[i]*e;
#else
    int max_n=n/8;
    DECLARE_REG_BI_COMPLEX(reg_c);
    REG_SPLAT_BI_COMPLEX(reg_c,c);
    DECLARE_REG_BI_COMPLEX(reg_e);
    REG_SPLAT_BI_COMPLEX(reg_e,e);
    NISSA_PARALLEL_LOOP(i,0,max_n)
      {
	DECLARE_REG_BI_HALFSPIN(reg_in1);
	DECLARE_REG_BI_HALFSPIN(reg_in2);
	DECLARE_REG_BI_HALFSPIN(reg_out);
	
	double *in1=(double*)(((bi_halfspin*)b)[i]);
	double *in2=(double*)(((bi_halfspin*)d)[i]);
	
	BI_HALFSPIN_PREFETCH_NEXT_NEXT(in1);
	BI_HALFSPIN_PREFETCH_NEXT_NEXT(in2);
	
	REG_LOAD_BI_HALFSPIN(reg_in1,in1);
	REG_LOAD_BI_HALFSPIN(reg_in2,in2);
	REG_BI_HALFSPIN_PROD_4DOUBLE(reg_out,reg_in1,reg_c);
	REG_BI_HALFSPIN_SUMM_THE_PROD_4DOUBLE(reg_out,reg_out,reg_in2,reg_e);
	STORE_REG_BI_HALFSPIN(((bi_halfspin*)a)[i],reg_out);
      }
    //last part
    if(IS_MASTER_THREAD) for(int i=max_n*8;i<n;i++) a[i]=b[i]*c+d[i]*e;
#endif  
    if(!(OPT&DO_NOT_SET_FLAGS)) set_borders_invalid(a);
  }}

  //////////////////////////////////////////////////////// quadruple precision ///////////////////////////////////////////
  
  //a=b
  THREADABLE_FUNCTION_3ARG(double_vector_from_quadruple_vector, double*,a, float_128*,b, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) a[i]=double_from_float_128(b[i]);set_borders_invalid(a);}}
  THREADABLE_FUNCTION_3ARG(quadruple_vector_from_double_vector, float_128*,a, double*,b, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) float_128_from_64(a[i],b[i]);set_borders_invalid(a);}}
  
  //a=a+b
  THREADABLE_FUNCTION_3ARG(quadruple_vector_summassign_double_vector, float_128*,a, double*,b, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) float_128_summassign_64(a[i],b[i]);set_borders_invalid(a);}}
  
  //a=b-c
  THREADABLE_FUNCTION_4ARG(quadruple_vector_subt_from_double_vector, float_128*,a, double*,b, float_128*,c, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) float_128_subt_from_64(a[i],b[i],c[i]);set_borders_invalid(a);}}
  
  /////////////////// scalar prodcut in quadruple /////////////////

  //(a,b)
  THREADABLE_FUNCTION_4ARG(quadruple_vector_glb_scalar_prod, float_128*,glb_res, float_128*,a, float_128*,b, int,n)
  {
    //perform thread summ
    float_128 loc_thread_res={0,0};
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(i,0,n)
      float_128_summ_the_prod(loc_thread_res,a[i],b[i]);
    
    glb_reduce_float_128(*glb_res,loc_thread_res);
  }}

  //(a,b)
  void double_conv_quadruple_vector_glb_scalar_prod(double *out,float_128 *a,float_128 *b,int n)
  {
    float_128 out_128;
    quadruple_vector_glb_scalar_prod(&out_128,a,b,n);
    *out=out_128[0];
  }

  //////////////// only quadruple accumulation //////////////////
  
  //(a,b)
  THREADABLE_FUNCTION_4ARG(quadruple_accumulate_double_vector_glb_scalar_prod, float_128*,a, double*,b, double*,c ,int,n)
  {
    //perform thread summ
    float_128 loc_thread_res={0,0};
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(i,0,n)
      float_128_summassign_64(loc_thread_res,b[i]*c[i]);
    
    glb_reduce_float_128(*a,loc_thread_res);
  }}

  //(a,b)
  double double_conv_quadruple_accumulate_double_vector_glb_scalar_prod(double *a,double *b,int n)
  {float_128 out;quadruple_accumulate_double_vector_glb_scalar_prod(&out,a,b,n);return double_from_float_128(out);}
  
  //////////////// color put/get from colorspinspin////////////////////////

  THREADABLE_FUNCTION_4ARG(get_color_from_colorspinspin, color*,out, colorspinspin*,in, int,id1, int,id2)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) get_color_from_colorspinspin(out[ivol],in[ivol],id1,id2);
    set_borders_invalid(out);
  }}

  THREADABLE_FUNCTION_4ARG(put_color_into_colorspinspin, colorspinspin*,out, color*,in, int,id1, int,id2)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) put_color_into_colorspinspin(out[ivol],in[ivol],id1,id2);
    set_borders_invalid(out);
  }}

  //////////////// color put/get from spincolor//////////////////

  THREADABLE_FUNCTION_3ARG(get_color_from_spincolor, color*,out, spincolor*,in, int,id)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) get_color_from_spincolor(out[ivol],in[ivol],id);
    set_borders_invalid(out);
  }}

  THREADABLE_FUNCTION_3ARG(put_color_into_spincolor, spincolor*,out, color*,in, int,id)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) put_color_into_spincolor(out[ivol],in[ivol],id);
    set_borders_invalid(out);
  }}

  //////////////// colorspinspin put/get ////////////////////////

  THREADABLE_FUNCTION_3ARG(get_spincolor_from_colorspinspin, spincolor*,out, colorspinspin*,in, int,id)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) get_spincolor_from_colorspinspin(out[ivol],in[ivol],id);
    set_borders_invalid(out);
  }}

  THREADABLE_FUNCTION_3ARG(put_spincolor_into_colorspinspin, colorspinspin*,out, spincolor*,in, int,id)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) put_spincolor_into_colorspinspin(out[ivol],in[ivol],id);
    set_borders_invalid(out);
  }}

  //////////////// su3spinspin put/get ////////////////////////
  
  THREADABLE_FUNCTION_4ARG(get_spincolor_from_su3spinspin, spincolor*,out, su3spinspin*,in, int,id, int,ic)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) get_spincolor_from_su3spinspin(out[ivol],in[ivol],id,ic);
    set_borders_invalid(out);
  }}

  THREADABLE_FUNCTION_4ARG(put_spincolor_into_su3spinspin, su3spinspin*,out, spincolor*,in, int,id, int,ic)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) put_spincolor_into_su3spinspin(out[ivol],in[ivol],id,ic);
    set_borders_invalid(out);
  }}

  ////////////////// spincolor algebra/////////////////////

  THREADABLE_FUNCTION_3ARG(safe_dirac_prod_spincolor, spincolor*,out, dirac_matr*,m, spincolor*,in)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) safe_dirac_prod_spincolor(out[ivol],m,in[ivol]);
    set_borders_invalid(out);
  }}

  ///////////////////// rotations ////////////////////////

  THREADABLE_FUNCTION_3ARG(rotate_vol_colorspinspin_to_physical_basis, colorspinspin*,s, int,rsi, int,rso)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int ic=0;ic<3;ic++)
	rotate_spinspin_to_physical_basis(s[ivol][ic],rsi,rso);
    set_borders_invalid(s);
  }}

  THREADABLE_FUNCTION_3ARG(rotate_vol_su3spinspin_to_physical_basis, su3spinspin*,s, int,rsi, int,rso)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int ic1=0;ic1<3;ic1++)
	for(int ic2=0;ic2<3;ic2++)
	  rotate_spinspin_to_physical_basis(s[ivol][ic1][ic2],rsi,rso);
    set_borders_invalid(s);
  }}

  THREADABLE_FUNCTION_3ARG(parallel_memcpy,void*,out, void*,in, int,n)
  {
#ifdef USE_THREADS
    GET_THREAD_ID();
    
    NISSA_CHUNK_WORKLOAD(start,chunk_load,end,0,n,thread_id,NACTIVE_THREADS);
    memcpy((char*)out+start,(char*)in+start,chunk_load);
    end++;//to avoid warning
    THREAD_BARRIER();
#else
    memcpy(out,in,n);
#endif
  }}
}
