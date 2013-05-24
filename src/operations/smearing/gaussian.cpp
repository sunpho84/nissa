#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/global_variables.h"
#include "../../communicate/communicate.h"
#include "../../base/vectors.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../routines/ios.h"
#include "../../routines/thread.h"

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
