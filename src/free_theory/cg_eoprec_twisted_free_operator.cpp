#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"
#include "free_theory_types.hpp"
#include "twisted_free_Dirac_eoprec_operator.hpp"

namespace nissa
{
  //invert Koo defined in equation (7)
  void inv_tmDkern_eoprec_square_eos(spin *sol,spin *guess,tm_quark_info qu,int nitermax,double residue,spin *source)
  {
    
    int niter=nitermax;
    int riter=0;
    int rniter=5;
    spin *p=nissa_malloc("p",locVolh+bord_volh,spin);
    spin *r=nissa_malloc("r",locVolh,spin);
    spin *s=nissa_malloc("s",locVolh,spin);
    spin *temp1=nissa_malloc("temp1",locVolh+bord_volh,spin);
    spin *temp2=nissa_malloc("temp2",locVolh+bord_volh,spin);
    
    ///////////////// prepare the internal source /////////////////
    
    if(guess==NULL) vector_reset(sol);
    else vector_copy(sol,guess);
    
    const int n=locVolh*sizeof(spin)/sizeof(double);
    
    //external loop, used if the internal exceed the maximal number of iterations
    double lambda; //(r_(k+1),r_(k+1))
    double source_norm;
    double_vector_glb_scalar_prod(&source_norm,(double*)source,(double*)source,n);
    do
      {
	//calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
	double delta;
	tmDkern_eoprec_square_eos(s,temp1,temp2,qu,sol);
	
	double_vector_subt((double*)r,(double*)source,(double*)s,n);
	double_vector_copy((double*)p,(double*)r,n);
	double_vector_glb_scalar_prod(&delta,(double*)r,(double*)r,n);
	
	if(riter==0)
	  {
	    master_printf("\nSource norm: %lg\n",source_norm);
	    master_printf("iter 0 relative residue: %lg\n",delta/source_norm);
	  }
	
	//main loop
	int iter=0;
	do
	  {
	    double omega; //(r_k,r_k)/(p_k*DD*p_k)
	    double alpha;
	    
	    tmDkern_eoprec_square_eos(s,temp1,temp2,qu,p);
	    double_vector_glb_scalar_prod(&alpha,(double*)s,(double*)p,n);
	    omega=delta/alpha;
	    
	    double_vector_summ_double_vector_prod_double((double*)sol,(double*)sol,(double*)p,omega,n);
	    double_vector_summ_double_vector_prod_double((double*)r,(double*)r,(double*)s,-omega,n);
	    double_vector_glb_scalar_prod(&lambda,(double*)r,(double*)r,n);
	    
	    double gammag=lambda/delta;
	    delta=lambda;
	    
	    //p_(k+1)=r_(k+1)+gammag*p_k
	    double_vector_summ_double_vector_prod_double((double*)p,(double*)r,(double*)p,gammag,n);
	    
	    iter++;
	    
	    if(iter%10==0) master_printf("iter %d relative residue: %lg\n",iter,lambda/source_norm);
	  }
	while(lambda>(residue*source_norm) && iter<niter);
	
	//last calculation of residual, in the case iter>niter
	tmDkern_eoprec_square_eos(s,temp1,temp2,qu,sol);
	double_vector_subt((double*)r,(double*)source,(double*)s,n);
	double_vector_glb_scalar_prod(&lambda,(double*)r,(double*)r,n);
	
	master_printf("\nfinal relative residue (after %d iters): %lg where %lg was required\n",iter,lambda/source_norm,residue);
	
	riter++;
      }
    while(lambda>(residue*source_norm) && riter<rniter);
    
    nissa_free(s);
    nissa_free(p);
    nissa_free(r);
    nissa_free(temp1);
    nissa_free(temp2);
  }
  
  //Invert twisted mass operator using e/o preconditioning.
  void inv_tmD_cg_eoprec_eos(spin *solution_lx,spin *guess_Koo,tm_quark_info qu,int nitermax,double residue,spin *source_lx)
  {
    
    //prepare the e/o split version of the source
    eo_ptr<spin> source_eos;
    source_eos[0]=nissa_malloc("source_eos0",locVolh+bord_volh,spin);
    source_eos[1]=nissa_malloc("source_eos1",locVolh+bord_volh,spin);
    split_lx_vector_into_eo_parts(source_eos,source_lx);
    
    //prepare the e/o split version of the solution
    eo_ptr<spin> solution_eos;
    solution_eos[0]=nissa_malloc("solution_eos_0",locVolh+bord_volh,spin);
    solution_eos[1]=nissa_malloc("solution_eos_0",locVolh+bord_volh,spin);
    
    ///////////////////////////////////// invert with e/o improvement ///////////////////////////////////
    
    spin *varphi=nissa_malloc("varphi",locVolh+bord_volh,spin);
    
    //Equation (8.a)
    spin *temp=nissa_malloc("temp",locVolh+bord_volh,spin);
    inv_tmDee_or_oo_eos(temp,qu,source_eos[EVN]);
    
    //Equation (8.b)
    tmn2Doe_eos(varphi,temp,qu.bc);
    NISSA_PARALLEL_LOOP(ieo,0,locVolh)
      for(int id=0;id<2;id++)
	for(int ri=0;ri<2;ri++)
	  { //gamma5 is explicitely wrote
	    varphi[ieo][id  ][ri]=+source_eos[ODD][ieo][id  ][ri]+varphi[ieo][id  ][ri]*0.5;
	    varphi[ieo][id+2][ri]=-source_eos[ODD][ieo][id+2][ri]-varphi[ieo][id+2][ri]*0.5;
	  }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(varphi);
    
    //Equation (9) using solution_eos[EVN] as temporary vector
    inv_tmDkern_eoprec_square_eos(temp,guess_Koo,qu,nitermax,residue,varphi);
    tm_quark_info mqu=qu;
    mqu.mass*=-1;
    tmDkern_eoprec_eos(solution_eos[ODD],solution_eos[EVN],mqu,temp);
    if(guess_Koo!=NULL) vector_copy(guess_Koo,temp); //if a guess was passed, return new one
    nissa_free(temp);
    
    //Equation (10)
    tmn2Deo_eos(varphi,solution_eos[ODD],qu.bc);
    NISSA_PARALLEL_LOOP(ieo,0,locVolh)
      for(int id=0;id<4;id++)
	for(int ri=0;ri<2;ri++)
	  varphi[ieo][id][ri]=source_eos[EVN][ieo][id][ri]+varphi[ieo][id][ri]*0.5;
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(varphi);
    inv_tmDee_or_oo_eos(solution_eos[EVN],qu,varphi);
    
    nissa_free(varphi);
    
    /////////////////////////// paste the e/o parts of the solution together and free ///////////////////
    
    paste_eo_parts_into_lx_vector(solution_lx,solution_eos);
    
    for(int par=0;par<2;par++)
      {
	nissa_free(source_eos[par]);
	nissa_free(solution_eos[par]);
      }
  }
}
