// Template to invert using bi-cg-stab
// macro to be defined:

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/bench.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#if CG_NARG >= 6
 #error not supported
#endif

namespace nissa
{

  //move elsewhere
  int nbicgstab_inv=0;
  double bicgstab_inv_over_time=0;
  

#if CG_NARG == 0
  THREADABLE_FUNCTION_5ARG(BICGSTAB_INVERT, BASETYPE*,sol, BASETYPE*,guess, int,niter, double,residue, BASETYPE*,source)
#elif CG_NARG == 1
  THREADABLE_FUNCTION_6ARG(BICGSTAB_INVERT, BASETYPE*,sol, BASETYPE*,guess, AT1,A1, int,niter, double,residue, BASETYPE*,source)
#elif CG_NARG == 2
  THREADABLE_FUNCTION_7ARG(BICGSTAB_INVERT, BASETYPE*,sol, BASETYPE*,guess, AT1,A1, AT2,A2, int,niter, double,residue, BASETYPE*,source)
#elif CG_NARG == 3
  THREADABLE_FUNCTION_8ARG(BICGSTAB_INVERT, BASETYPE*,sol, BASETYPE*,guess, AT1,A1, AT2,A2, AT3,A3, int,niter, double,residue, BASETYPE*,source)
#elif CG_NARG == 4
  THREADABLE_FUNCTION_9ARG(BICGSTAB_INVERT, BASETYPE*,sol, BASETYPE*,guess, AT1,A1, AT2,A2, AT3,A3, AT4,A4, int,niter, double,residue, BASETYPE*,source)
#elif CG_NARG == 5
  THREADABLE_FUNCTION_10ARG(BICGSTAB_INVERT, BASETYPE*,sol, BASETYPE*,guess, AT1,A1, AT2,A2, AT3,A3, AT4,A4, AT5,A5, int,niter, double,residue, BASETYPE*,source)
#endif
  {
    GET_THREAD_ID();

    verbosity_lv2_master_printf("\n");

    int riter=0;
    BASETYPE *t=nissa_malloc("s",BULK_VOL,BASETYPE);
    BASETYPE *s=nissa_malloc("s",BULK_VOL,BASETYPE);
    BASETYPE *p=nissa_malloc("p",BULK_VOL+BORD_VOL,BASETYPE);
    BASETYPE *nu=nissa_malloc("nu",BULK_VOL,BASETYPE);
    BASETYPE *r=nissa_malloc("r",BULK_VOL,BASETYPE);
    BASETYPE *rhat0=nissa_malloc("rhat0",BULK_VOL,BASETYPE);
    
    //macro to be defined externally, allocating all the required additional vectors
    CG_ADDITIONAL_VECTORS_ALLOCATION();
    if(guess==NULL) vector_reset(sol);
    else vector_copy(sol,guess);
    
    START_TIMING(bicgstab_inv_over_time,nbicgstab_inv);
    
    int each=VERBOSITY_LV3?1:10;
    
    //reset p and nu
    vector_reset(p);
    vector_reset(nu);
    
    //compute the norm of the source
    double source_norm;
    double_vector_glb_scalar_prod(&source_norm,(double*)source,(double*)source,BULK_VOL*NDOUBLES_PER_SITE);

    //calculate rhat0=r0=b-DD*sol_0 and delta_0=(p0,p0)
    APPLY_OPERATOR(s,CG_OPERATOR_PARAMETERS sol);    
    double_vector_subt((double*)r,(double*)source,(double*)s,BULK_VOL*NDOUBLES_PER_SITE);
    vector_copy(rhat0,r);
    double delta;
    double_vector_glb_scalar_prod(&delta,(double*)r,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
    
    //print source info and check it
    if(riter==0) verbosity_lv2_master_printf("Source norm: %lg\n",source_norm);
    if(source_norm==0 || std::isnan(source_norm)) crash("invalid norm: %lg",source_norm);
    verbosity_lv2_master_printf("iter 0 relative residue: %lg\n",delta/source_norm);
    
    int final_iter;
    
    //main loop
    int iter=0;
    double alpha=1,omega=1,lambda,rho=1;
    do
      {
	//this is already iter 1
	final_iter=(++iter);
	
	double rho_old=rho,omega_old=omega;
	  
	//compute (rhat0,r) and beta
	double_vector_glb_scalar_prod(&rho,(double*)rhat0,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
	double beta=(rho/rho_old)*(alpha/omega_old);
	
	//compute p=r+beta*(p-omega_old*nu)
	double_vector_summ_double_vector_prod_double((double*)p,(double*)p,(double*)nu,-omega_old,BULK_VOL*NDOUBLES_PER_SITE);
	double_vector_summ_double_vector_prod_double((double*)p,(double*)r,(double*)p,beta,BULK_VOL*NDOUBLES_PER_SITE);
	
	//nu=A*P
	if(IS_MASTER_THREAD) bicgstab_inv_over_time+=take_time();
	APPLY_OPERATOR(nu,CG_OPERATOR_PARAMETERS p);
	if(IS_MASTER_THREAD) bicgstab_inv_over_time-=take_time();
	
	//compute alpha=rho/(rhat0,nu)
	double temp;
	double_vector_glb_scalar_prod(&temp,(double*)rhat0,(double*)nu,BULK_VOL*NDOUBLES_PER_SITE);
	alpha=rho/temp;
	
	//s=r-alpha*nu
	double_vector_summ_double_vector_prod_double((double*)s,(double*)r,(double*)nu,-alpha,BULK_VOL*NDOUBLES_PER_SITE);
	//t=As

	if(IS_MASTER_THREAD) bicgstab_inv_over_time+=take_time();
	APPLY_OPERATOR(t,CG_OPERATOR_PARAMETERS s);
	if(IS_MASTER_THREAD) bicgstab_inv_over_time-=take_time();
	
	//omega=(t,s)/(t,t)
	double_vector_glb_scalar_prod(&omega,(double*)t,(double*)s,BULK_VOL*NDOUBLES_PER_SITE);
	double_vector_glb_scalar_prod(&temp,(double*)t,(double*)t,BULK_VOL*NDOUBLES_PER_SITE);
	omega/=temp;
	
	//x=x+alpha*p+omega*s
	double_vector_summ_double_vector_prod_double((double*)sol,(double*)sol,(double*)p,alpha,BULK_VOL*NDOUBLES_PER_SITE);
	double_vector_summ_double_vector_prod_double((double*)sol,(double*)sol,(double*)s,omega,BULK_VOL*NDOUBLES_PER_SITE);
	
	//r=s-t*omega	
	double_vector_summ_double_vector_prod_double((double*)r,(double*)s,(double*)t,-omega,BULK_VOL*NDOUBLES_PER_SITE);
	
	//lambda=(r,r)
	double_vector_glb_scalar_prod(&lambda,(double*)r,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
	
	if(iter%each==0) verbosity_lv2_master_printf("iter %d relative residue: %lg\n",iter,lambda/source_norm);
      }
    while(lambda>=(residue*source_norm) && iter<niter);
    
    //last calculation of residual
    APPLY_OPERATOR(s,CG_OPERATOR_PARAMETERS sol);
    double_vector_subt((double*)r,(double*)source,(double*)s,BULK_VOL*NDOUBLES_PER_SITE);
    double_vector_glb_scalar_prod(&lambda,(double*)r,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
    
    verbosity_lv2_master_printf("final relative residue (after %d iters): %lg where %lg was required\n",
				final_iter,lambda/source_norm,residue);

    //check if not converged
    if(final_iter==niter) crash("exit without converging");
    
    if(IS_MASTER_THREAD) bicgstab_inv_over_time+=take_time();
    
    nissa_free(t);
    nissa_free(s);
    nissa_free(p);
    nissa_free(nu);
    nissa_free(r);
    nissa_free(rhat0);
    
    //macro to be defined externally
    CG_ADDITIONAL_VECTORS_FREE();
  }
  THREADABLE_FUNCTION_END
}
