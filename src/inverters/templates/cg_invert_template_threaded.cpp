// Template to invert using c.g.
// macro to be defined:
//  -apply_operator
//  -cg_invert
//  -cg_parameters_proto
//  -cg_inner_parameters_call
//  -BULK_VOL, BORD_VOL
//  -basetype
//  -NDOUBLES_PER_SITE

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/bench.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

#if CG_NARG >= 6
 #error not supported
#endif

namespace nissa
{
#if CG_NARG == 0
  void CG_INVERT(BASETYPE* sol,BASETYPE* guess,int niter,double residue,BASETYPE* source)
#elif CG_NARG == 1
  void CG_INVERT(BASETYPE* sol,BASETYPE* guess,AT1 A1,int niter,double residue,BASETYPE* source)
#elif CG_NARG == 2
  void CG_INVERT(BASETYPE* sol,BASETYPE* guess,AT1 A1,AT2 A2,int niter,double residue,BASETYPE* source)
#elif CG_NARG == 3
  void CG_INVERT(BASETYPE* sol,BASETYPE* guess,AT1 A1,AT2 A2,AT3 A3,int niter,double residue,BASETYPE* source)
#elif CG_NARG == 4
  void CG_INVERT(BASETYPE* sol,BASETYPE* guess,AT1 A1,AT2 A2,AT3 A3,AT4 A4,int niter,double residue,BASETYPE* source)
#elif CG_NARG == 5
  void CG_INVERT(BASETYPE* sol,BASETYPE* guess,AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,int niter,double residue,BASETYPE* source)
#endif
  {
    
    verbosity_lv2_master_printf("\n");
    
    BASETYPE *s=nissa_malloc("s",BULK_VOL,BASETYPE);
    BASETYPE *p=nissa_malloc("p",BULK_VOL+BORD_VOL,BASETYPE);
    BASETYPE *r=nissa_malloc("r",BULK_VOL,BASETYPE);
    
    //macro to be defined externally, allocating all the required additional vectors
    CG_ADDITIONAL_VECTORS_ALLOCATION();
    if(guess==NULL) vector_reset(sol);
    else vector_copy(sol,guess);
    
    START_TIMING(cg_inv_over_time,ncg_inv);
    int each=VERBOSITY_LV3?1:10;
    
    double source_norm;
    double_vector_glb_scalar_prod(&source_norm,(double*)source,(double*)source,BULK_VOL*NDOUBLES_PER_SITE);
    //calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
    APPLY_OPERATOR(s,CG_OPERATOR_PARAMETERS sol);
    
    double_vector_subt((double*)r,(double*)source,(double*)s,BULK_VOL*NDOUBLES_PER_SITE);
    double_vector_copy((double*)p,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
    double delta;
    double_vector_glb_scalar_prod(&delta,(double*)r,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
    
    verbosity_lv2_master_printf("Source norm: %lg\n",source_norm);
    if(source_norm==0 || std::isnan(source_norm)) crash("invalid norm: %lg",source_norm);
    verbosity_lv2_master_printf("iter 0 relative residue: %lg\n",delta/source_norm);
    
    int final_iter;
    
    //main loop
    int iter=0;
    double alpha,omega,gammag,lambda;
    do
      {
	//this is already iter 1
	final_iter=(++iter);
	
	//(r_k,r_k)/(p_k*DD*p_k)
	STOP_TIMING(cg_inv_over_time);
	APPLY_OPERATOR(s,CG_OPERATOR_PARAMETERS p);
	cg_inv_over_time-=take_time();
	
	double_vector_glb_scalar_prod(&alpha,(double*)s,(double*)p,BULK_VOL*NDOUBLES_PER_SITE);
	omega=delta/alpha;
	
	//sol_(k+1)=x_k+omega*p_k
	double_vector_summ_double_vector_prod_double((double*)sol,(double*)sol,(double*)p,omega,BULK_VOL*NDOUBLES_PER_SITE);
	//r_(k+1)=x_k-omega*p_k
	double_vector_summ_double_vector_prod_double((double*)r,(double*)r,(double*)s,-omega,BULK_VOL*NDOUBLES_PER_SITE);
	//(r_(k+1),r_(k+1))
	double_vector_glb_scalar_prod(&lambda,(double*)r,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
	
	//(r_(k+1),r_(k+1))/(r_k,r_k)
	gammag=lambda/delta;
	delta=lambda;
	
	//checks
	if(std::isnan(gammag)) crash("nanned");
	
	//p_(k+1)=r_(k+1)+gammag*p_k
	double_vector_summ_double_vector_prod_double((double*)p,(double*)r,(double*)p,gammag,BULK_VOL*NDOUBLES_PER_SITE);
	
	if(iter%each==0) verbosity_lv2_master_printf("iter %d relative residue: %lg\n",iter,lambda/source_norm);
      }
    while(lambda>=(residue*source_norm) && iter<niter);
    
    //last calculation of residual
    APPLY_OPERATOR(s,CG_OPERATOR_PARAMETERS sol);
    double_vector_subt((double*)r,(double*)source,(double*)s,BULK_VOL*NDOUBLES_PER_SITE);
    double_vector_glb_scalar_prod(&lambda,(double*)r,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
    
    verbosity_lv2_master_printf("final relative residue (after %d iters): %lg where %lg was required\n",
				final_iter,lambda/source_norm,residue);
    if(lambda/source_norm>=2*residue)
      master_printf("WARNING: true residue %lg much larger than required and expected one %lg\n",
		    lambda/source_norm,residue);
    
    verbosity_lv1_master_printf(" Total cg iterations: %d\n",final_iter);
    
    //check if not converged
    if(final_iter==niter) crash("exit without converging");
    
    cg_inv_over_time+=take_time();
    
    nissa_free(s);
    nissa_free(p);
    nissa_free(r);
    
    //macro to be defined externally
    CG_ADDITIONAL_VECTORS_FREE();
  }
}

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
