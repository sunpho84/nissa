#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <math.h>

#include "communicate/communicate.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "new_types/complex.hpp"
#include "new_types/dirac.hpp"
#include "new_types/float_128.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"
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
    NISSA_PARALLEL_LOOP(i,0,n) a[i]=0;
    
    set_borders_invalid(a);
  }
  
  //copy
  template <class T1,class T2> void internal_vector_copy(T1 *a,T2 *b,int n)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(i,0,n) a[i]=b[i];
    
    set_borders_invalid(a);
  }
  
  //double to single and vv
  THREADABLE_FUNCTION_3ARG(double_vector_to_single, float*,a, double*,b, int,n)
  {internal_vector_copy(a,b,n);}THREADABLE_FUNCTION_END
  THREADABLE_FUNCTION_3ARG(single_vector_to_double, double*,a, float*,b, int,n)
  {internal_vector_copy(a,b,n);}THREADABLE_FUNCTION_END
  //double to double
  THREADABLE_FUNCTION_3ARG(double_vector_copy, double*,a, double*,b, int,n)
  {internal_vector_copy(a,b,n);}THREADABLE_FUNCTION_END
  
  //set to zero
  void single_vector_init_to_zero(float *a,int n)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(i,0,n) a[i]=0;
    
    set_borders_invalid(a);
  }
  
  //copy
  void single_vector_copy(float *a,float *b,int n)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(i,0,n) a[i]=b[i];
    
    set_borders_invalid(a);
  }
  
  //summ
  THREADABLE_FUNCTION_4ARG(double_vector_summ, double*,out, double*,in1, double*,in2, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) out[i]=in1[i]+in2[i];set_borders_invalid(out);}THREADABLE_FUNCTION_END
  
  //subt
  THREADABLE_FUNCTION_4ARG(double_vector_subt, double*,out, double*,in1, double*,in2, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) out[i]=in1[i]-in2[i];set_borders_invalid(out);}THREADABLE_FUNCTION_END
  
  //prod with double
  THREADABLE_FUNCTION_4ARG(double_vector_prod_double, double*,out, double*,in, double,r, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) out[i]=r*in[i];set_borders_invalid(out);}THREADABLE_FUNCTION_END
  THREADABLE_FUNCTION_4ARG(float_128_vector_prod_double, float_128*,out, float_128*,in, double,r, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) float_128_prod_64(out[i],in[i],r);set_borders_invalid(out);}THREADABLE_FUNCTION_END
  
  //prod with double of the summ
  THREADABLE_FUNCTION_5ARG(double_vector_prod_the_summ_double, double*,out, double,r, double*,in1, double*,in2, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) out[i]=r*(in1[i]+in2[i]);set_borders_invalid(out);}THREADABLE_FUNCTION_END
  
  //scalar product
  THREADABLE_FUNCTION_4ARG(double_vector_glb_scalar_prod, double*,glb_res, double*,a, double*,b, int,n)
  {
#ifndef REPRODUCIBLE_RUN
    //perform thread summ
    double loc_thread_res=0;
    GET_THREAD_ID();
    
#ifndef BGQ
    NISSA_PARALLEL_LOOP(i,0,n) loc_thread_res+=a[i]*b[i];
#else
    int max_n=n/8;
    DECLARE_REG_VIR_HALFSPIN(reg_loc_thread_res);
    REG_SPLAT_VIR_HALFSPIN(reg_loc_thread_res,0);
    NISSA_PARALLEL_LOOP(i,0,max_n)
      {
	DECLARE_REG_VIR_HALFSPIN(reg_a);
	DECLARE_REG_VIR_HALFSPIN(reg_b);
	
	double *a_ptr=(double*)(((vir_halfspin*)a)[i]);
	double *b_ptr=(double*)(((vir_halfspin*)b)[i]);
	
	VIR_HALFSPIN_PREFETCH_NEXT_NEXT(a_ptr);
	VIR_HALFSPIN_PREFETCH_NEXT_NEXT(b_ptr);
	REG_LOAD_VIR_HALFSPIN(reg_a,a_ptr);
	REG_LOAD_VIR_HALFSPIN(reg_b,b_ptr);
	REG_VIR_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(reg_loc_thread_res,s0),NAME2(reg_loc_thread_res,s0),NAME2(reg_a,s0),NAME2(reg_b,s0));
	REG_VIR_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(reg_loc_thread_res,s1),NAME2(reg_loc_thread_res,s1),NAME2(reg_a,s1),NAME2(reg_b,s1));
      }
    REG_VIR_COMPLEX_SUMM(reg_loc_thread_res_s0,reg_loc_thread_res_s0,reg_loc_thread_res_s1);
    
    //store
    vir_complex temp;
    STORE_REG_VIR_COMPLEX(temp,reg_loc_thread_res_s0);
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
  }
  THREADABLE_FUNCTION_END
  
  //scalar product
  THREADABLE_FUNCTION_4ARG(single_vector_glb_scalar_prod, float*,glb_res, float*,a, float*,b, int,n)
  {
    //perform thread summ
    float loc_thread_res=0;
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(i,0,n)
      loc_thread_res+=a[i]*b[i];
    
    (*glb_res)=glb_reduce_single(loc_thread_res);
  }
  THREADABLE_FUNCTION_END
  
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
  }
  THREADABLE_FUNCTION_END
  
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
  }
  THREADABLE_FUNCTION_END
  
  //put a vector of complex equal to its conjugate
  THREADABLE_FUNCTION_3ARG(complex_vector_conj, complex*,res, complex*,in, int,n)
  {
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(i,0,n) complex_conj(res[i],in[i]);
    set_borders_invalid(res);
  }
  THREADABLE_FUNCTION_END
  
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
  }
  THREADABLE_FUNCTION_END
  
  //a[]=b[]+c[]*d
  THREADABLE_FUNCTION_6ARG(double_vector_summ_double_vector_prod_double, double*,a, double*,b, double*,c, double,d, int,n, int,OPT)
  {
    GET_THREAD_ID();
#ifndef BGQ
    NISSA_PARALLEL_LOOP(i,0,n) a[i]=b[i]+c[i]*d;
#else
    int max_n=n/8;
    DECLARE_REG_VIR_COMPLEX(reg_d);
    REG_SPLAT_VIR_COMPLEX(reg_d,d);
    NISSA_PARALLEL_LOOP(i,0,max_n)
      {
	DECLARE_REG_VIR_HALFSPIN(reg_in1);
	DECLARE_REG_VIR_HALFSPIN(reg_in2);
	DECLARE_REG_VIR_HALFSPIN(reg_out);
	
	double *in1=(double*)(((vir_halfspin*)b)[i]);
	double *in2=(double*)(((vir_halfspin*)c)[i]);
	
	VIR_HALFSPIN_PREFETCH_NEXT_NEXT(in1);
	VIR_HALFSPIN_PREFETCH_NEXT_NEXT(in2);
	
	REG_LOAD_VIR_HALFSPIN(reg_in1,in1);
	REG_LOAD_VIR_HALFSPIN(reg_in2,in2);
	REG_VIR_HALFSPIN_SUMM_THE_PROD_4DOUBLE(reg_out,reg_in1,reg_in2,reg_d);
	STORE_REG_VIR_HALFSPIN(((vir_halfspin*)a)[i],reg_out);
      }
    //last part
    if(IS_MASTER_THREAD) for(int i=max_n*8;i<n;i++) a[i]=b[i]+c[i]*d;
#endif
    if(!(OPT&DO_NOT_SET_FLAGS)) set_borders_invalid(a);
  }
  THREADABLE_FUNCTION_END
  
  //single version
  THREADABLE_FUNCTION_6ARG(single_vector_summ_single_vector_prod_single, float*,a, float*,b, float*,c, float,d, int,n, int,OPT)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(i,0,n) a[i]=b[i]+c[i]*d;
    
    if(!(OPT&DO_NOT_SET_FLAGS)) set_borders_invalid(a);
  }
  THREADABLE_FUNCTION_END

  //a[]=b[]*c+d[]*e
  THREADABLE_FUNCTION_7ARG(double_vector_linear_comb, double*,a, double*,b, double,c, double*,d, double,e, int,n, int,OPT)
  {
    GET_THREAD_ID();
#ifndef BGQ
    NISSA_PARALLEL_LOOP(i,0,n) a[i]=b[i]*c+d[i]*e;
#else
    int max_n=n/8;
    DECLARE_REG_VIR_COMPLEX(reg_c);
    REG_SPLAT_VIR_COMPLEX(reg_c,c);
    DECLARE_REG_VIR_COMPLEX(reg_e);
    REG_SPLAT_VIR_COMPLEX(reg_e,e);
    NISSA_PARALLEL_LOOP(i,0,max_n)
      {
	DECLARE_REG_VIR_HALFSPIN(reg_in1);
	DECLARE_REG_VIR_HALFSPIN(reg_in2);
	DECLARE_REG_VIR_HALFSPIN(reg_out);
	
	double *in1=(double*)(((vir_halfspin*)b)[i]);
	double *in2=(double*)(((vir_halfspin*)d)[i]);
	
	VIR_HALFSPIN_PREFETCH_NEXT_NEXT(in1);
	VIR_HALFSPIN_PREFETCH_NEXT_NEXT(in2);
	
	REG_LOAD_VIR_HALFSPIN(reg_in1,in1);
	REG_LOAD_VIR_HALFSPIN(reg_in2,in2);
	REG_VIR_HALFSPIN_PROD_4DOUBLE(reg_out,reg_in1,reg_c);
	REG_VIR_HALFSPIN_SUMM_THE_PROD_4DOUBLE(reg_out,reg_out,reg_in2,reg_e);
	STORE_REG_VIR_HALFSPIN(((vir_halfspin*)a)[i],reg_out);
      }
    //last part
    if(IS_MASTER_THREAD) for(int i=max_n*8;i<n;i++) a[i]=b[i]*c+d[i]*e;
#endif  
    if(!(OPT&DO_NOT_SET_FLAGS)) set_borders_invalid(a);
  }
  THREADABLE_FUNCTION_END

  //a[]=b[]*c+d[]*e
  THREADABLE_FUNCTION_7ARG(single_vector_linear_comb, float*,a, float*,b, float,c, float*,d, float,e, int,n, int,OPT)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(i,0,n) a[i]=b[i]*c+d[i]*e;
    
    if(!(OPT&DO_NOT_SET_FLAGS)) set_borders_invalid(a);
  }
  THREADABLE_FUNCTION_END

  ////////////////////////////////////////////////////// quadruple precision ///////////////////////////////////////////
  
  //a=b
  THREADABLE_FUNCTION_3ARG(double_vector_from_quadruple_vector, double*,a, float_128*,b, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) a[i]=double_from_float_128(b[i]);set_borders_invalid(a);}
  THREADABLE_FUNCTION_END
  THREADABLE_FUNCTION_3ARG(quadruple_vector_from_double_vector, float_128*,a, double*,b, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) float_128_from_64(a[i],b[i]);set_borders_invalid(a);}
  THREADABLE_FUNCTION_END
  
  //a=a+b
  THREADABLE_FUNCTION_3ARG(quadruple_vector_summassign_double_vector, float_128*,a, double*,b, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) float_128_summassign_64(a[i],b[i]);set_borders_invalid(a);}
  THREADABLE_FUNCTION_END
  
  //a=b-c
  THREADABLE_FUNCTION_4ARG(quadruple_vector_subt_from_double_vector, float_128*,a, double*,b, float_128*,c, int,n)
  {GET_THREAD_ID();NISSA_PARALLEL_LOOP(i,0,n) float_128_subt_from_64(a[i],b[i],c[i]);set_borders_invalid(a);}
  THREADABLE_FUNCTION_END
  
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
  }
  THREADABLE_FUNCTION_END

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
  }
  THREADABLE_FUNCTION_END
  
  //(a,b)
  double double_conv_quadruple_accumulate_double_vector_glb_scalar_prod(double *a,double *b,int n)
  {float_128 out;quadruple_accumulate_double_vector_glb_scalar_prod(&out,a,b,n);return double_from_float_128(out);}
  
  ////////////////////// color put/get from su3 ///////////////////////////
  
  THREADABLE_FUNCTION_3ARG(get_color_from_su3, color**,out, su3**,in, int,ic_source)
  {
    GET_THREAD_ID();
    for(int eo=0;eo<2;eo++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh) get_color_from_su3(out[eo][ivol],in[eo][ivol],ic_source);
	set_borders_invalid(out[eo]);
      }
  }
  THREADABLE_FUNCTION_END
  
    THREADABLE_FUNCTION_3ARG(put_color_into_su3, su3**,out, color**,in, int,ic_source)
  {
    GET_THREAD_ID();
    for(int eo=0;eo<2;eo++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh) put_color_into_su3(out[eo][ivol],in[eo][ivol],ic_source);
	set_borders_invalid(out[eo]);
      }
  }
  THREADABLE_FUNCTION_END
  
  //////////////// color put/get from colorspinspin////////////////////////
  
  THREADABLE_FUNCTION_4ARG(get_color_from_colorspinspin, color*,out, colorspinspin*,in, int,id1, int,id2)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) get_color_from_colorspinspin(out[ivol],in[ivol],id1,id2);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_4ARG(put_color_into_colorspinspin, colorspinspin*,out, color*,in, int,id1, int,id2)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) put_color_into_colorspinspin(out[ivol],in[ivol],id1,id2);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //////////////// color put/get from spincolor//////////////////
  
  THREADABLE_FUNCTION_3ARG(get_color_from_spincolor, color*,out, spincolor*,in, int,id)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) get_color_from_spincolor(out[ivol],in[ivol],id);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END

  THREADABLE_FUNCTION_3ARG(put_color_into_spincolor, spincolor*,out, color*,in, int,id)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) put_color_into_spincolor(out[ivol],in[ivol],id);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //////////////// colorspinspin put/get ////////////////////////
  
  THREADABLE_FUNCTION_3ARG(get_spincolor_from_colorspinspin, spincolor*,out, colorspinspin*,in, int,id)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) get_spincolor_from_colorspinspin(out[ivol],in[ivol],id);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_3ARG(put_spincolor_into_colorspinspin, colorspinspin*,out, spincolor*,in, int,id)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) put_spincolor_into_colorspinspin(out[ivol],in[ivol],id);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //////////////// su3spinspin put/get ////////////////////////
  
  THREADABLE_FUNCTION_4ARG(get_spincolor_from_su3spinspin, spincolor*,out, su3spinspin*,in, int,id, int,ic)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) get_spincolor_from_su3spinspin(out[ivol],in[ivol],id,ic);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_4ARG(put_spincolor_into_su3spinspin, su3spinspin*,out, spincolor*,in, int,id, int,ic)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) put_spincolor_into_su3spinspin(out[ivol],in[ivol],id,ic);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  ////////////////// spincolor algebra/////////////////////
  
  THREADABLE_FUNCTION_3ARG(safe_dirac_prod_spincolor, spincolor*,out, dirac_matr*,m, spincolor*,in)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) safe_dirac_prod_spincolor(out[ivol],m,in[ivol]);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END

  THREADABLE_FUNCTION_3ARG(safe_dirac_prod_colorspinspin, colorspinspin*,out, dirac_matr*,m, colorspinspin*,in)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) safe_dirac_prod_colorspinspin(out[ivol],m,in[ivol]);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  ///////////////////// rotations ////////////////////////
  
  THREADABLE_FUNCTION_3ARG(rotate_vol_colorspinspin_to_physical_basis, colorspinspin*,s, int,rsi, int,rso)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int ic=0;ic<3;ic++)
	rotate_spinspin_to_physical_basis(s[ivol][ic],rsi,rso);
    set_borders_invalid(s);
  }
  THREADABLE_FUNCTION_END

  THREADABLE_FUNCTION_3ARG(rotate_vol_su3spinspin_to_physical_basis, su3spinspin*,s, int,rsi, int,rso)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int ic1=0;ic1<3;ic1++)
	for(int ic2=0;ic2<3;ic2++)
	  rotate_spinspin_to_physical_basis(s[ivol][ic1][ic2],rsi,rso);
    set_borders_invalid(s);
  }
  THREADABLE_FUNCTION_END
  
  //ildg to nissa and vice-versa
  THREADABLE_FUNCTION_1ARG(quad_su3_nissa_to_ildg_reord_in_place, quad_su3*,in)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      quad_su3_nissa_to_ildg_reord(in[ivol],in[ivol]);
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
  THREADABLE_FUNCTION_1ARG(quad_su3_ildg_to_nissa_reord_in_place, quad_su3*,in)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      quad_su3_ildg_to_nissa_reord(in[ivol],in[ivol]);
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END

  
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
  }
  THREADABLE_FUNCTION_END
}
