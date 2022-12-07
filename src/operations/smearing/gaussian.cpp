#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/field.hpp>
#include <base/vectors.hpp>
#include <communicate/borders.hpp>
#include <geometry/geometry_lx.hpp>
#include <linalgs/linalgs.hpp>
#include <new_types/su3_op.hpp>
#include <routines/ios.hpp>

namespace nissa
{
  //apply kappa*H
#define DEFINE_GAUSSIAN_SMEARING_APPLY_KAPPA_H(TYPE)			\
  void NAME2(gaussian_smearing_apply_kappa_H,TYPE)(LxField<TYPE>& H,	\
						   const double* kappa,	\
						   const LxField<quad_su3>& conf, \
						   const LxField<TYPE>& in) \
  {									\
    in.updateHalo();							\
    conf.updateHalo();							\
    									\
    H.reset();								\
									\
    NISSA_PARALLEL_LOOP(ivol,0,locVol)					\
      {									\
	for(int mu=1;mu<NDIM;mu++)					\
	  {								\
	    const int ivup=loclxNeighup[ivol][mu];			\
	    const int ivdw=loclxNeighdw[ivol][mu];			\
	    TYPE temp;							\
									\
	    NAME2(unsafe_su3_prod,TYPE)(temp,conf[ivol][mu],in[ivup]);	\
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
  void gaussian_smearing(LxField<TYPE>& smear_sc,			\
			 const LxField<TYPE>& origi_sc,			\
			 const LxField<quad_su3>& conf,			\
			 const double* kappa,				\
			 const int& niter,				\
			 LxField<TYPE>* ext_temp,			\
			 LxField<TYPE>* ext_H)				\
  {									\
    if(niter<1)								\
      {									\
	verbosity_lv2_master_printf("Skipping smearing (0 iter required)\n"); \
	if(smear_sc!=origi_sc) smear_sc=origi_sc;			\
      }									\
    else								\
      {									\
	LxField<TYPE> *_temp=ext_temp;					\
	if(_temp==nullptr) _temp=new LxField<TYPE>("temp",WITH_HALO);	\
	LxField<TYPE>& temp=*_temp;					\
									\
	LxField<TYPE> *_H=ext_H;					\
	if(ext_H==nullptr) _H=new LxField<TYPE>("H",WITH_HALO);	\
	LxField<TYPE>& H=*_H;						\
  									\
	const double norm_fact=1/(1+2*(kappa[1]+kappa[2]+kappa[3]));	\
									\
	verbosity_lv2_master_printf("GAUSSIAN smearing with kappa={%g,%g,%g}, %d iterations\n",kappa[1],kappa[2],kappa[3],niter); \
									\
	/*iter 0*/							\
	temp=origi_sc;							\
									\
	/*loop over gaussian iterations*/				\
	for(int iter=0;iter<niter;iter++)				\
	  {								\
	    verbosity_lv3_master_printf("GAUSSIAN smearing with kappa={%g,%g,%g} iteration %d of %d\n",kappa[1],kappa[2],kappa[3],iter,niter); \
									\
	    /*apply kappa*H*/						\
	    NAME2(gaussian_smearing_apply_kappa_H,TYPE)(H,kappa,conf,temp); \
	    /*add kappa*H and dynamic normalize*/			\
	    temp+=H;							\
	    temp*=norm_fact;						\
	  }								\
									\
	smear_sc=temp;							\
									\
	if(ext_H==nullptr) delete _H;					\
	if(ext_temp==nullptr) delete _temp;				\
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
  
#undef DEFINE_GAUSSIAN_SMEARING
#undef DEFINE_GAUSSIAN_SMEARING_APPLY_KAPPA_H
  
  //smear with a polynomial of H
  template <typename TYPE>
  void gaussian_smearing(LxField<TYPE>& smear_sc,
			 const LxField<TYPE>& origi_sc,
			 const LxField<quad_su3>& conf,
			 const double& kappa,
			 const int& nterm,
			 const double *coeff,
			 const int *exponent)
  {
    if(nterm==0 or (nterm==1 and exponent[0]==0 and coeff[0]==1))
      {
	if(smear_sc!=origi_sc) smear_sc==origi_sc;
      }
    else
      {
	//copy to a temp buffer
	LxField<TYPE> temp1("temp1",WITH_HALO);
	temp1=origi_sc;
	
	//allocate two temps vectors for gaussian
	LxField<TYPE> temp2("temp2",WITH_HALO);
	LxField<TYPE> temp3("temp3",WITH_HALO);
	
	//reset the output
	smear_sc.reset();
	
	for(int iterm=0;iterm<nterm;iterm++)
	  {
	    //compute the number of smearing steps
	    int nstep=exponent[iterm];
	    if(iterm>0)
	      nstep-=exponent[iterm-1];
	    
	    //smear
	    gaussian_smearing(temp1,temp1,conf,kappa,nstep,temp2,temp3);
	    
	    //accumulate
	    temp1*=coeff[iterm];
	    smear_sc+=temp1;
	    //improvable
	  }
	
	set_borders_invalid(smear_sc);
      }
  }
}
