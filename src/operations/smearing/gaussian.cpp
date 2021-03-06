#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //apply kappa*H
#define DEFINE_GAUSSIAN_SMEARING_APPLY_KAPPA_H(TYPE)			\
  void NAME2(gaussian_smearing_apply_kappa_H,TYPE)(TYPE* H,double* kappa,quad_su3* conf,TYPE* in) \
  {									\
									\
    NAME3(communicate_lx,TYPE,borders)(in);				\
    communicate_lx_quad_su3_borders(conf);				\
									\
    vector_reset(H);							\
									\
    NISSA_PARALLEL_LOOP(ivol,0,locVol)					\
      {									\
	for(int mu=1;mu<NDIM;mu++)					\
	  {								\
	    int ivup=loclxNeighup[ivol][mu];				\
	    int ivdw=loclxNeighdw[ivol][mu];				\
	    TYPE temp;							\
	    NAME2(unsafe_su3_prod,TYPE)(temp,conf[ivol][mu],in[ivup]); \
	    NAME2(su3_dag_summ_the_prod,TYPE)(temp,conf[ivdw][mu],in[ivdw]);	\
	    NAME2(TYPE,summ_the_prod_double)(H[ivol],temp,kappa[mu]);		\
	  }								\
      }									\
    NISSA_PARALLEL_LOOP_END;						\
    									\
    set_borders_invalid(H);						\
  }									\

//gaussian smearing
#define DEFINE_GAUSSIAN_SMEARING(TYPE)					\
  void gaussian_smearing(TYPE* smear_sc,TYPE* origi_sc,quad_su3* conf,double* kappa,int niter,TYPE* ext_temp,TYPE* ext_H) \
  {									\
    if(niter<1)								\
      {									\
	verbosity_lv2_master_printf("Skipping smearing (0 iter required)\n"); \
	if(smear_sc!=origi_sc) vector_copy(smear_sc,origi_sc);		\
      }									\
    else								\
      {									\
	TYPE *temp=ext_temp;						\
	if(temp==NULL) temp=nissa_malloc("temp",locVol+bord_vol,TYPE); \
									\
	TYPE *H=ext_H;							\
	if(ext_H==NULL) H=nissa_malloc("H",locVol+bord_vol,TYPE);	\
									\
	double norm_fact=1/(1+2*(kappa[1]+kappa[2]+kappa[3]));		\
									\
	verbosity_lv2_master_printf("GAUSSIAN smearing with kappa={%g,%g,%g}, %d iterations\n",kappa[1],kappa[2],kappa[3],niter); \
									\
	/*iter 0*/							\
	double_vector_copy((double*)temp,(double*)origi_sc,locVol*sizeof(TYPE)/sizeof(double)); \
									\
	/*loop over gaussian iterations*/				\
	for(int iter=0;iter<niter;iter++)				\
	  {								\
	    verbosity_lv3_master_printf("GAUSSIAN smearing with kappa={%g,%g,%g} iteration %d of %d\n",kappa[1],kappa[2],kappa[3],iter,niter); \
									\
	    /*apply kappa*H*/						\
	    NAME2(gaussian_smearing_apply_kappa_H,TYPE)(H,kappa,conf,temp); \
	    /*add kappa*H and dynamic normalize*/			\
	    double_vector_prod_the_summ_double((double*)temp,norm_fact,(double*)temp,(double*)H,sizeof(TYPE)/sizeof(double)*locVol); \
	  }								\
									\
	double_vector_copy((double*)smear_sc,(double*)temp,locVol*sizeof(TYPE)/sizeof(double)); \
									\
	if(ext_H==NULL) nissa_free(H);					\
	if(ext_temp==NULL) nissa_free(temp);				\
      }									\
  }									\
  
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
	TYPE *temp1=nissa_malloc("temp",locVol+bord_vol,TYPE);
	vector_copy(temp1,origi_sc);
	
	//allocate two temp vectors for gaussian
	TYPE *temp2=nissa_malloc("temp2",locVol+bord_vol,TYPE);
	TYPE *temp3=nissa_malloc("temp3",locVol+bord_vol,TYPE);
	
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
	    double_vector_summ_double_vector_prod_double((double*)smear_sc,(double*)smear_sc,(double*)temp1,coeff[iterm],locVol*sizeof(TYPE)/sizeof(double));
	  }
	
	set_borders_invalid(smear_sc);
	
	nissa_free(temp1);
	nissa_free(temp2);
	nissa_free(temp3);
      }
  }
}
