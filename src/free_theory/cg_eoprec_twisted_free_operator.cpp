#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif
#include "free_theory_types.hpp"
#include "twisted_free_Dirac_eoprec_operator.hpp"

namespace nissa
{
  //invert Koo defined in equation (7)
  void inv_tmDkern_eoprec_square_eos(spin *sol,spin *guess,tm_quark_info qu,int nitermax,double residue,spin *source)
  {
    GET_THREAD_ID();
    
    int niter=nitermax;
    int riter=0;
    int rniter=5;
    spin *p=nissa_malloc("p",loc_volh+bord_volh,spin);
    spin *r=nissa_malloc("r",loc_volh,spin);
    spin *s=nissa_malloc("s",loc_volh,spin);
    spin *temp1=nissa_malloc("temp1",loc_volh+bord_volh,spin);
    spin *temp2=nissa_malloc("temp2",loc_volh+bord_volh,spin);
    
    ///////////////// prepare the internal source /////////////////
    
    if(guess==NULL) vector_reset(sol);
    else vector_copy(sol,guess);
    
    //external loop, used if the internal exceed the maximal number of iterations
    double lambda; //(r_(k+1),r_(k+1))
    double source_norm;
    do
      {
	//calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
	double delta;
	{
	  tmDkern_eoprec_square_eos(s,temp1,temp2,qu,sol);
	  
	  double loc_delta=0,loc_source_norm=0;
	  NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	    for(int id=0;id<4;id++)
	      for(int ri=0;ri<2;ri++)
		{
		  double c1=source[ivol][id][ri]-s[ivol][id][ri];
		  p[ivol][id][ri]=r[ivol][id][ri]=c1;
		  if(riter==0) loc_source_norm+=source[ivol][id][ri]*source[ivol][id][ri];
		  loc_delta+=c1*c1;
		}
	  NISSA_PARALLEL_LOOP_END;
	  set_borders_invalid(p);
	  delta=glb_reduce_double(loc_delta);
	  
	  if(riter==0)
	    {
	      source_norm=glb_reduce_double(loc_source_norm);
	      master_printf("\nSource norm: %lg\n",source_norm);
	      master_printf("iter 0 relative residue: %lg\n",delta/source_norm);
	    }
	}
	
	//main loop
	int iter=0;
	do
	  {
	    double omega; //(r_k,r_k)/(p_k*DD*p_k)
	    double alpha;
	    
	    tmDkern_eoprec_square_eos(s,temp1,temp2,qu,p);
	    
	    double loc_alpha=0; //real part of the scalar product
	    NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	      for(int id=0;id<4;id++)
		for(int ri=0;ri<2;ri++)
		  loc_alpha+=s[ivol][id][ri]*p[ivol][id][ri];
	    NISSA_PARALLEL_LOOP_END;
	    alpha=glb_reduce_double(loc_alpha);
	    omega=delta/alpha;
	    
	    double loc_lambda=0;
	    NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	      for(int id=0;id<4;id++)
		for(int ri=0;ri<2;ri++)
		  {
		    sol[ivol][id][ri]+=omega*p[ivol][id][ri];    //sol_(k+1)=x_k+omega*p_k
		    double c1=r[ivol][id][ri]-omega*s[ivol][id][ri];//r_(k+1)=x_k-omega*pk
		    r[ivol][id][ri]=c1;
		    loc_lambda+=c1*c1;
		  }
	    NISSA_PARALLEL_LOOP_END;
	    set_borders_invalid(sol);
	    lambda=glb_reduce_double(loc_lambda);
	    
	    double gammag=lambda/delta; //(r_(k+1),r_(k+1))/(r_k,r_k)
	    delta=lambda;
	    
	    //p_(k+1)=r_(k+1)+gammag*p_k
	    NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	      for(int id=0;id<4;id++)
		for(int ri=0;ri<2;ri++)
		  p[ivol][id][ri]=r[ivol][id][ri]+gammag*p[ivol][id][ri];
	    NISSA_PARALLEL_LOOP_END;
	    set_borders_invalid(p);
	    
	    iter++;
	    
	    if(iter%10==0) master_printf("iter %d relative residue: %lg\n",iter,lambda/source_norm);
	  }
	while(lambda>(residue*source_norm) && iter<niter);
	
	//last calculation of residual, in the case iter>niter
	tmDkern_eoprec_square_eos(s,temp1,temp2,qu,sol);
	{
	  double loc_lambda=0;
	  NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	    for(int id=0;id<4;id++)
	      for(int ri=0;ri<2;ri++)
		{
		  double c1=source[ivol][id][ri]-s[ivol][id][ri];
		  loc_lambda+=c1*c1;
		}
	  NISSA_PARALLEL_LOOP_END;
	  
	  lambda=glb_reduce_double(loc_lambda);
	}
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
    GET_THREAD_ID();
    
    //prepare the e/o split version of the source
    spin *source_eos[2];
    source_eos[0]=nissa_malloc("source_eos0",loc_volh+bord_volh,spin);
    source_eos[1]=nissa_malloc("source_eos1",loc_volh+bord_volh,spin);
    split_lx_vector_into_eo_parts(source_eos,source_lx);
    
    //prepare the e/o split version of the solution
    spin *solution_eos[2];
    solution_eos[0]=nissa_malloc("solution_eos_0",loc_volh+bord_volh,spin);
    solution_eos[1]=nissa_malloc("solution_eos_0",loc_volh+bord_volh,spin);
    
    ///////////////////////////////////// invert with e/o improvement ///////////////////////////////////
    
    spin *varphi=nissa_malloc("varphi",loc_volh+bord_volh,spin);
    
    //Equation (8.a)
    spin *temp=nissa_malloc("temp",loc_volh+bord_volh,spin);
    inv_tmDee_or_oo_eos(temp,qu,source_eos[EVN]);
    
    //Equation (8.b)
    tmn2Doe_eos(varphi,temp,qu.bc);
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
      for(int id=0;id<2;id++)
	for(int ri=0;ri<2;ri++)
	  { //gamma5 is explicitely wrote
	    varphi[ivol][id  ][ri]=+source_eos[ODD][ivol][id  ][ri]+varphi[ivol][id  ][ri]*0.5;
	    varphi[ivol][id+2][ri]=-source_eos[ODD][ivol][id+2][ri]-varphi[ivol][id+2][ri]*0.5;
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
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
      for(int id=0;id<4;id++)
	for(int ri=0;ri<2;ri++)
	  varphi[ivol][id][ri]=source_eos[EVN][ivol][id][ri]+varphi[ivol][id][ri]*0.5;
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
