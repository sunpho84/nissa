#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_Wsklx.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //apply kappa*H
#define DEFINE_GAUSSIAN_SMEARING_APPLY_KAPPA_H(TYPE)			\
  THREADABLE_FUNCTION_4ARG(NAME2(gaussian_smearing_apply_kappa_H,TYPE), TYPE*,H, double,kappa, quad_su3*,conf, TYPE*,in) \
  {									\
    GET_THREAD_ID();							\
									\
    NAME3(communicate_lx,TYPE,borders)(in);				\
    communicate_lx_quad_su3_borders(conf);				\
									\
    vector_reset(H);							\
									\
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)					\
      {									\
	for(int mu=1;mu<4;mu++)						\
	  {								\
	    int ivup=loclx_neighup[ivol][mu];				\
	    int ivdw=loclx_neighdw[ivol][mu];				\
									\
	    NAME2(su3_summ_the_prod,TYPE)(H[ivol],conf[ivol][mu],in[ivup]); \
	    NAME2(su3_dag_summ_the_prod,TYPE)(H[ivol],conf[ivdw][mu],in[ivdw]); \
	  }								\
	NAME2(TYPE,prod_double)(H[ivol],H[ivol],kappa);			\
      }									\
									\
    set_borders_invalid(H);						\
  }}

//gaussian smearing
#define DEFINE_GAUSSIAN_SMEARING(TYPE)					\
  THREADABLE_FUNCTION_7ARG(gaussian_smearing, TYPE*,smear_sc, TYPE*,origi_sc, quad_su3*,conf, double,kappa, int,niter, TYPE*,ext_temp, TYPE*,ext_H) \
  {									\
    if(niter<1)								\
      {									\
	verbosity_lv2_master_printf("Skipping smearing (0 iter required)\n"); \
	if(smear_sc!=origi_sc) vector_copy(smear_sc,origi_sc);		\
      }									\
    else								\
      {									\
	TYPE *temp=ext_temp;						\
	if(temp==NULL) temp=nissa_malloc("temp",loc_vol+bord_vol,TYPE); \
									\
	TYPE *H=ext_H;							\
	if(ext_H==NULL) H=nissa_malloc("H",loc_vol+bord_vol,TYPE);	\
									\
	double norm_fact=1/(1+6*kappa);					\
									\
	verbosity_lv2_master_printf("GAUSSIAN smearing with kappa=%g, %d iterations\n",kappa,niter); \
									\
	/*iter 0*/							\
	vector_copy(temp,origi_sc);					\
									\
	/*loop over gaussian iterations*/				\
	for(int iter=0;iter<niter;iter++)				\
	  {								\
	    verbosity_lv3_master_printf("GAUSSIAN smearing with kappa=%g iteration %d of %d\n",kappa,iter,niter); \
									\
	    /*apply kappa*H*/						\
	    NAME2(gaussian_smearing_apply_kappa_H,TYPE)(H,kappa,conf,temp); \
	    /*add kappa*H and dynamic normalize*/			\
	    double_vector_prod_the_summ_double((double*)temp,norm_fact,(double*)temp,(double*)H,sizeof(TYPE)/sizeof(double)*loc_vol); \
	  }								\
									\
	vector_copy(smear_sc,temp);					\
									\
	if(ext_H==NULL) nissa_free(H);					\
	if(ext_temp==NULL) nissa_free(temp);				\
      }									\
  }}

DEFINE_GAUSSIAN_SMEARING_APPLY_KAPPA_H(su3spinspin)
DEFINE_GAUSSIAN_SMEARING(su3spinspin)

DEFINE_GAUSSIAN_SMEARING_APPLY_KAPPA_H(colorspinspin)
DEFINE_GAUSSIAN_SMEARING(colorspinspin)

DEFINE_GAUSSIAN_SMEARING_APPLY_KAPPA_H(spincolor)
DEFINE_GAUSSIAN_SMEARING(spincolor)

DEFINE_GAUSSIAN_SMEARING_APPLY_KAPPA_H(color)
DEFINE_GAUSSIAN_SMEARING(color)

//smear with a polynomial of H
template <class TYPE> void gaussian_smearing(TYPE *smear_sc,TYPE *origi_sc,quad_su3 *conf,double kappa,int nterm,double *coeff,int *exponent)
{
  if(nterm==0||(nterm==1&&exponent[0]==0&&coeff[0]==1)){if(smear_sc!=origi_sc) vector_copy(smear_sc,origi_sc);}
  else
    {
      //copy to a temp buffer
      TYPE *temp1=nissa_malloc("temp",loc_vol+bord_vol,TYPE);
      vector_copy(temp1,origi_sc);
      
      //allocate two temp vectors for gaussian
      TYPE *temp2=nissa_malloc("temp2",loc_vol+bord_vol,TYPE);
      TYPE *temp3=nissa_malloc("temp3",loc_vol+bord_vol,TYPE);
      
      //reset the output
      vector_reset(smear_sc);
      
      for(int iterm=0;iterm<nterm;iterm++)
	{
	  //compute the number of smearing steps
	  int nstep=exponent[iterm];
	  if(iterm>0) nstep-=exponent[iterm-1];
	  
	  //smear
	  gaussian_smearing(temp1,temp1,conf,kappa,nstep,temp2,temp3);
	  
	  //accumulate
	  double_vector_summ_double_vector_prod_double((double*)smear_sc,(double*)smear_sc,(double*)temp1,coeff[iterm],loc_vol*sizeof(TYPE)/sizeof(double));
	}
      
      set_borders_invalid(smear_sc);
      
      nissa_free(temp1);
      nissa_free(temp2);
      nissa_free(temp3);
    }
}

//apply just the eight point terms
#define DEFINE_APPLY_DER(T)						\
  void apply_der(T *temp,oct_su3 *conf,T *in,int ivol)			\
  {									\
    int *inter_fr_in_pos=Wsklx_hopping_matrix_output_pos.inter_fr_in_pos+ivol*8; \
    for(int mu=0;mu<4;mu++)						\
      NAME2(unsafe_su3_prod,T)(temp[inter_fr_in_pos[0+mu]],conf[ivol][mu],in[ivol]); \
    for(int mu=0;mu<4;mu++)						\
      NAME2(unsafe_su3_dag_prod,T)(temp[inter_fr_in_pos[4+mu]],conf[ivol][4+mu],in[ivol]); \
  }

DEFINE_APPLY_DER(spincolor)
DEFINE_APPLY_DER(colorspinspin)
DEFINE_APPLY_DER(su3spinspin)

template <class T> void gaussian_smearing_iter(T *out,T *in,oct_su3 *conf,double kappa,T *temp,comm_t &comm)
{
  GET_THREAD_ID();
  
  if(get_vect(temp)->nel<8*loc_vol+bord_vol) crash("allocate border for temp vector");
  
#ifndef BGQCACCA
  const int ndouble_per_site=sizeof(T)/sizeof(double);
#else
  int nbi_complex_per_site=sizeof(T)/sizeof(bi_complex);
  if(sizeof(T)%sizeof(bi_complex)) crash("not good for BGQ!");
#endif
  
  //apply derivative on surface and blast it
  NISSA_PARALLEL_LOOP(isurf,0,surf_vol) apply_der(temp,conf,in,isurf);
  THREAD_BARRIER();
  parallel_memcpy(send_buf,temp+8*loc_vol,bord_vol*sizeof(T));
  comm_start(comm);
  
  //apply derivative on bulk and wait
  NISSA_PARALLEL_LOOP(ibulk,surf_vol,loc_vol) apply_der(temp,conf,in,ibulk);
  comm_wait(comm);
  
  //put incoming pieces to temp
  NISSA_PARALLEL_LOOP(ipiece,0,surf_vol)
    memcpy(temp[Wsklx_hopping_matrix_output_pos.final_fr_inter_pos[ipiece]],((T*)recv_buf)[ipiece],sizeof(T));
  THREAD_BARRIER();
  
  //summ the eight piecese together with original vector and put normalization
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    {
      //#ifndef BGQ
      double *a=(double*)(out[ivol]),*b=(double*)(temp[ivol*8+0]),*c=(double*)(in[ivol]);
      double acc[ndouble_per_site];
      //summ 8 of the 6 pieces: 0 and 4 must be jumped
      for(int i=0;i<ndouble_per_site;i++) acc[i]=b[ndouble_per_site*1+i]+b[ndouble_per_site*2+i];
      for(int i=0;i<ndouble_per_site;i++) acc[i]+=b[ndouble_per_site*3+i];
      for(int ip=5;ip<8;ip++) for(int i=0;i<ndouble_per_site;i++) acc[i]+=b[ndouble_per_site*ip+i];
      //normalize
      double norm_fact=1/(1+6*kappa);
      for(int i=0;i<ndouble_per_site;i++) a[i]=norm_fact*(c[i]+kappa*acc[i]);
    }
  //#else
  //#endif
  
  set_borders_invalid(out);
}

#define DEFINE_GAUSSIAN_SMEARING_SINK_BASED(T,COMM)			\
  THREADABLE_FUNCTION_6ARG(gaussian_smearing, T*,out, T*,in, oct_su3*,conf, double,kappa, T*,temp, int,niter) \
  {									\
if(niter==0) {if(out!=in) vector_copy(out,in);}				\
 else									\
   {									\
     gaussian_smearing_iter(out,in,conf,kappa,temp,COMM);		\
     for(int iter=1;iter<niter;iter++)					\
       gaussian_smearing_iter(out,out,conf,kappa,temp,COMM);		\
   }									\
}}									\
void gaussian_smearing_sink_based(T *ext_out,T *ext_in,quad_su3 *ext_conf,double kappa,int niter,T *aux_temp=NULL,T *aux_out=NULL,T *aux_in=NULL,oct_su3 *aux_conf=NULL) \
{									\
  T *out=(aux_out!=NULL)?aux_out:nissa_malloc("out",loc_vol,T);		\
  T *temp=(aux_temp!=NULL)?aux_temp:nissa_malloc("temp",loc_vol*8+bord_vol,T); \
  T *in=(aux_in!=NULL)?aux_in:nissa_malloc("in",loc_vol,T);		\
  oct_su3 *conf=(aux_conf!=NULL)?aux_conf:nissa_malloc("conf",loc_vol,oct_su3);	\
  									\
  lx_remap_to_Wsklx(in,ext_in);						\
  lx_conf_remap_to_Wsklx(conf,ext_conf);				\
  gaussian_smearing(out,in,conf,kappa,temp,niter);			\
  Wsklx_remap_to_lx(ext_out,out);					\
  									\
  if(aux_conf==NULL) nissa_free(conf);					\
  if(aux_in==NULL) nissa_free(in);					\
  if(aux_temp==NULL) nissa_free(temp);					\
  if(aux_out==NULL) nissa_free(out);					\
}

DEFINE_GAUSSIAN_SMEARING_SINK_BASED(spincolor,lx_spincolor_comm)
DEFINE_GAUSSIAN_SMEARING_SINK_BASED(colorspinspin,lx_colorspinspin_comm)
DEFINE_GAUSSIAN_SMEARING_SINK_BASED(su3spinspin,lx_su3spinspin_comm)
}
