// Template to invert using c.g.
// macro to be defined:
//  -apply_operator
//  -cg_invert
//  -cg_parameters_proto
//  -cg_inner_parameters_call
//  -BULK_VOL, BORD_VOL
//  -basetype
//  -NDOUBLES_PER_SITE

extern double cg_inv_over_time;
extern int ncg_inv;

#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../routines/ios.h"
#include "../../routines/openmp.h"

#if CGM_NARG >= 3
 #error not supported
#endif

#if CGM_NARG == 0
THREADABLE_FUNCTION_6ARG(CG_INVERT, BASETYPE*,sol, BASETYPE*,guess, int,niter, int,rniter, double,residue, BASETYPE*,source)
#elif CGM_NARG == 1
THREADABLE_FUNCTION_7ARG(CG_INVERT, BASETYPE*,sol, BASETYPE*,guess, AT1,A1, int,niter, int,rniter, double,residue, BASETYPE*,source)
#elif CGM_NARG == 2
THREADABLE_FUNCTION_8ARG(CG_INVERT, BASETYPE*,sol, BASETYPE*,guess, AT1,A1, AT2,A2, int,niter, int,rniter, double,residue, BASETYPE*,source)
#endif
{
  int riter=0;
  BASETYPE *s=nissa_malloc("s",BULK_VOL,BASETYPE);
  BASETYPE *p=nissa_malloc("p",BULK_VOL+BORD_VOL,BASETYPE);
  BASETYPE *r=nissa_malloc("r",BULK_VOL,BASETYPE);

  //macro to be defined externally, allocating all the required additional vectors
  CG_ADDITIONAL_VECTORS_ALLOCATION();
  if(guess==NULL) vector_reset(sol);
  else vector_copy(sol,guess);
  
  if(IS_MASTER_THREAD)
    {
      ncg_inv++;
      cg_inv_over_time-=take_time();
    }
  
  int each_list[4]={0,100,10,1},each;
  if(nissa_verbosity>=3) each=1;
  else each=each_list[nissa_verbosity];
  
  //external loop, used if the internal exceed the maximal number of iterations
  double source_norm,lambda;
  do
    {
      double_vector_glb_scalar_prod(&source_norm,(double*)source,(double*)source,BULK_VOL*NDOUBLES_PER_SITE);
      //calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
      APPLY_OPERATOR(s,CG_OPERATOR_PARAMETERS sol);
      
      double_vector_subt_double_vector_prod_double((double*)r,(double*)source,(double*)s,1,BULK_VOL*NDOUBLES_PER_SITE);
      double_vector_copy((double*)p,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
      double delta;
      double_vector_glb_scalar_prod(&delta,(double*)r,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
      
      if(riter==0) verbosity_lv2_master_printf("Source norm: %lg\n",source_norm);
      if(source_norm==0 || isnan(source_norm)) crash("invalid norm: %lg",source_norm);
      verbosity_lv2_master_printf("iter 0 relative residue: %lg\n",delta/source_norm);
      
      int final_iter;
      
      //main loop
      int iter=0;
      double alpha,omega,gammag,internal_lambda,internal_delta=delta;
      do
	{	  
	  //this is already iter 1
	  final_iter=(++iter);
	  
	  //(r_k,r_k)/(p_k*DD*p_k)
	  if(IS_MASTER_THREAD) cg_inv_over_time+=take_time();
	  APPLY_OPERATOR(s,CG_OPERATOR_PARAMETERS p);
	  if(IS_MASTER_THREAD) cg_inv_over_time-=take_time();
	  
	  double_vector_glb_scalar_prod(&alpha,(double*)s,(double*)p,BULK_VOL*NDOUBLES_PER_SITE);
	  omega=internal_delta/alpha;
	  
	  //sol_(k+1)=x_k+omega*p_k
	  double_vector_summ_double_vector_prod_double((double*)sol,(double*)sol,(double*)p,omega,BULK_VOL*NDOUBLES_PER_SITE);
	  //r_(k+1)=x_k-omega*p_k
	  double_vector_summ_double_vector_prod_double((double*)r,(double*)r,(double*)s,-omega,BULK_VOL*NDOUBLES_PER_SITE);
	  //(r_(k+1),r_(k+1))
	  double_vector_glb_scalar_prod(&internal_lambda,(double*)r,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
	  
	  //(r_(k+1),r_(k+1))/(r_k,r_k)
	  gammag=internal_lambda/internal_delta;
	  internal_delta=internal_lambda;
	  
	  //p_(k+1)=r_(k+1)+gammag*p_k
	  double_vector_summ_double_vector_prod_double((double*)p,(double*)r,(double*)p,gammag,BULK_VOL*NDOUBLES_PER_SITE);
	  
	  if(iter%each==0) verbosity_lv2_master_printf("iter %d relative residue: %lg\n",iter,internal_lambda/source_norm);
	}
      while(internal_lambda>(residue*source_norm) && iter<niter);
      
      //last calculation of residual, in the case iter>niter
      APPLY_OPERATOR(s,CG_OPERATOR_PARAMETERS sol);
      double_vector_subt_double_vector_prod_double((double*)r,(double*)source,(double*)s,1,BULK_VOL*NDOUBLES_PER_SITE);
      double_vector_glb_scalar_prod(&lambda,(double*)r,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
      
      verbosity_lv1_master_printf("\nfinal relative residue (after %d iters): %lg where %lg was required\n",final_iter,lambda/source_norm,residue);
      
      riter++;
    }
  while(lambda>(residue*source_norm) && riter<rniter);
  
  if(IS_MASTER_THREAD) cg_inv_over_time+=take_time();
  
  nissa_free(s);
  nissa_free(p);
  nissa_free(r);
  
  //macro to be defined externally
  CG_ADDITIONAL_VECTORS_FREE();
}}

#undef BASETYPE
#undef NDOUBLES_PER_SITE
#undef BULK_VOL
#undef BORD_VOL

#undef APPLY_OPERATOR
#undef CG_OPERATOR_PARAMETERS
#undef CG_INVERT
#undef CG_ADDITIONAL_VECTORS_FREE
#undef CG_ADDITIONAL_VECTORS_ALLOCATION
#undef CG_NARG
