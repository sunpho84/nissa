// Template to invert using c.g.
// macro to be defined:
//  -apply_operator
//  -cg_invert
//  -cg_parameters_proto
//  -cg_inner_parameters_call
//  -size_of_bulk, size_of_bord
//  -basetype
//  -ndoubles_per_site

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <omp.h>

#include "routines/ios.hpp"
namespace nissa
{
  void cg_invert(basetype *sol,basetype *guess,cg_parameters_proto,int niter,int rniter,double residue,basetype *source)
  {
    int riter=0;
    basetype *s=nissa_malloc("s",size_of_bulk,basetype);
    basetype *p=nissa_malloc("p",size_of_bulk+size_of_bord,basetype);
    basetype *r=nissa_malloc("r",size_of_bulk,basetype);
    
    //macro to be defined externally, allocating all the required additional vectors
    cg_additional_vectors_allocation();
    
    if(guess==NULL) vector_reset(sol);
    else vector_copy(sol,guess);
    
#ifdef BENCH
    ncg_inv++;
    cg_inv_over_time-=take_time();
#endif
    
    int each_list[4]={0,100,10,1},each;
    if(verbosity_lv>=3) each=1;
    else each=each_list[verbosity_lv];
    
    //external loop, used if the internal exceed the maximal number of iterations
    double source_norm,lambda;
    do
      {
	//calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
	apply_operator(s,cg_inner_parameters_call,sol);
	
	double_vector_subt((double*)r,(double*)source,(double*)s,size_of_bulk*ndoubles_per_site);
	double_vector_copy((double*)p,(double*)r,size_of_bulk*ndoubles_per_site);
	double_vector_glb_scalar_prod(&source_norm,(double*)source,(double*)source,size_of_bulk*ndoubles_per_site);
	double delta;
	double_vector_glb_scalar_prod(&delta,(double*)r,(double*)r,size_of_bulk*ndoubles_per_site);
	
	if(riter==0) verbosity_lv2_master_printf("Source norm: %lg\n",source_norm);
	if(source_norm==0 || isnan(source_norm)) crash("invalid norm: %lg",source_norm);
	
	verbosity_lv2_master_printf("iter 0 relative residue: %lg\n",delta/source_norm);
	
	int final_iter;
	
	//main loop
	int iter=0;
	double alpha,omega,gammag,internal_lambda,internal_delta=delta;
	do
	  {	  
	    //(r_k,r_k)/(p_k*DD*p_k)
#ifdef BENCH
	    cg_inv_over_time+=take_time();
#endif
	    apply_operator(s,cg_inner_parameters_call,p);
#ifdef BENCH
	    cg_inv_over_time-=take_time();
#endif
	    
	    double_vector_glb_scalar_prod(&alpha,(double*)s,(double*)p,size_of_bulk*ndoubles_per_site);
	    omega=internal_delta/alpha;
	    
	    //sol_(k+1)=x_k+omega*p_k
	    double_vector_summ_double_vector_prod_double((double*)sol,(double*)sol,(double*)p,omega,size_of_bulk*ndoubles_per_site);
	    //r_(k+1)=x_k-omega*p_k
	    double_vector_summ_double_vector_prod_double((double*)r,(double*)r,(double*)s,-omega,size_of_bulk*ndoubles_per_site);
	    //(r_(k+1),r_(k+1))
	    double_vector_glb_scalar_prod(&internal_lambda,(double*)r,(double*)r,size_of_bulk*ndoubles_per_site);
	    
	    //(r_(k+1),r_(k+1))/(r_k,r_k)
	    gammag=internal_lambda/internal_delta;
	    internal_delta=internal_lambda;
	    
	    //p_(k+1)=r_(k+1)+gammag*p_k
	    double_vector_summ_double_vector_prod_double((double*)p,(double*)r,(double*)p,gammag,size_of_bulk*ndoubles_per_site);
	    
	    final_iter=(++iter);
	    if(iter%each==0) verbosity_lv2_master_printf("iter %d relative residue: %lg\n",iter,internal_lambda/source_norm);
	  }
	while(internal_lambda>(residue*source_norm) && iter<niter);
	
	//last calculation of residual, in the case iter>niter
	apply_operator(s,cg_inner_parameters_call,sol);
	double_vector_subt((double*)r,(double*)source,(double*)s,size_of_bulk*ndoubles_per_site);
	double_vector_glb_scalar_prod(&lambda,(double*)r,(double*)r,size_of_bulk*ndoubles_per_site);
	
	verbosity_lv1_master_printf("\nfinal relative residue (after %d iters): %lg where %lg was required\n",final_iter,lambda/source_norm,residue);
	
	riter++;
      }
    while(lambda>(residue*source_norm) && riter<rniter);
    
#ifdef BENCH
    cg_inv_over_time+=take_time();
#endif
    
    nissa_free(s);
    nissa_free(p);
    nissa_free(r);
    
    //macro to be defined externally
    cg_additional_vectors_free();
  }
}

#undef basetype
#undef ndoubles_per_site
#undef size_of_bulk
#undef size_of_bord
  
#undef apply_operator
#undef cg_operator_parameters
#undef cg_invert
#undef cg_additional_vectors_free
#undef cg_additional_vectors_allocation
  
  
