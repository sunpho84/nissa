#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/bench.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  /*
    This is the prorotipe for a multi-shift inverter.
    The calls to the operator, the internal vectors definitions and the additional parameters must be defined thorugh macro.
    The file must be included inside another file defining all the macros.
    See "cgm_invert_tmQ2.c" as an example.
  */

#if CGM_NARG >= 6
#error not supported
#endif
  
#if CGM_NARG == 0
  THREADABLE_FUNCTION_6ARG(CGM_INVERT, BASETYPE**,sol, double*,shift, int,nshift, int,niter_max, double*,ext_req_res, BASETYPE*,source)
#elif CGM_NARG == 1
  THREADABLE_FUNCTION_7ARG(CGM_INVERT, BASETYPE**,sol, AT1,A1, double*,shift, int,nshift, int,niter_max, double*,ext_req_res, BASETYPE*,source)
#elif CGM_NARG == 2
  THREADABLE_FUNCTION_8ARG(CGM_INVERT, BASETYPE**,sol, AT1,A1, AT2,A2, double*,shift, int,nshift, int,niter_max, double*,ext_req_res, BASETYPE*,source)
#elif CGM_NARG == 3
  THREADABLE_FUNCTION_9ARG(CGM_INVERT, BASETYPE**,sol, AT1,A1, AT2,A2, AT3,A3, double*,shift, int,nshift, int,niter_max, double*,ext_req_res, BASETYPE*,source)
#elif CGM_NARG == 4
  THREADABLE_FUNCTION_10ARG(CGM_INVERT, BASETYPE**,sol, AT1,A1, AT2,A2, AT3,A3, AT4,A4, double*,shift, int,nshift, int,niter_max, double*,ext_req_res, BASETYPE*,source)
#elif CGM_NARG == 5
  THREADABLE_FUNCTION_11ARG(CGM_INVERT, BASETYPE**,sol, AT1,A1, AT2,A2, AT3,A3, AT4,A4, AT5,A5, double*,shift, int,nshift, int,niter_max, double*,ext_req_res, BASETYPE*,source)
#endif
  {
    GET_THREAD_ID();
    
#ifdef CG_128_INVERT
    //limit inner solver precision
    double max_inner_solver=1.0e-25;
    //used for inner solver in the case of 128 bit precision
    double inn_req_res[nshift];
    for(int ishift=0;ishift<nshift;ishift++)
      if(use_128_bit_precision && ext_req_res[ishift]<max_inner_solver)
	{
	  verbosity_lv2_master_printf("changing the inner solver residue for shift %d to %lg\n",ishift,max_inner_solver);
	  inn_req_res[ishift]=max_inner_solver;
	}
      else inn_req_res[ishift]=ext_req_res[ishift];
#else
    double *inn_req_res=ext_req_res;
#endif
    
#ifdef SQRT_SHIFT
    double sqrt_shift[nshift];
    for(int ishift=0;ishift<nshift;ishift++)
      sqrt_shift[ishift]=sqrt(shift[ishift]);
#define IN_SHIFT sqrt_shift
#else
#define IN_SHIFT shift
#endif
    
    if(IS_MASTER_THREAD)
      {
	ncgm_inv++;
	cgm_inv_over_time-=take_time();
      }
    
    int each=VERBOSITY_LV3?1:10;
    
    //macro to be defined externally, allocating all the required additional vectors
    CGM_ADDITIONAL_VECTORS_ALLOCATION();
    
    BASETYPE *s=nissa_malloc("s",BULK_VOL,BASETYPE);
    BASETYPE *r=nissa_malloc("r",BULK_VOL,BASETYPE);
    BASETYPE *p=nissa_malloc("p",BULK_VOL+BORD_VOL,BASETYPE);
    BASETYPE *ps[nshift];
    for(int ishift=0;ishift<nshift;ishift++) ps[ishift]=nissa_malloc("ps",BULK_VOL,BASETYPE);
    
    //     -sol[*]=0
    //     -ps[*]=source
    for(int ishift=0;ishift<nshift;ishift++)
      {
	double_vector_copy((double*)(ps[ishift]),(double*)source,BULK_VOL*NDOUBLES_PER_SITE);
	double_vector_init_to_zero((double*)(sol[ishift]),BULK_VOL*NDOUBLES_PER_SITE);
      }
    
    //     -p=source
    //     -r=source
    //     -calculate source_norm=(r,r)
    double_vector_copy((double*)p,(double*)source,BULK_VOL*NDOUBLES_PER_SITE);
    double_vector_copy((double*)r,(double*)source,BULK_VOL*NDOUBLES_PER_SITE);
    double source_norm;
    double_vector_glb_scalar_prod(&source_norm,(double*)r,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
    
    //writes source norm
    verbosity_lv2_master_printf(" Source norm: %lg\n",source_norm);
    if(source_norm==0 || std::isnan(source_norm)) crash("invalid norm: %lg",source_norm);
    
    //writes initial residue
    verbosity_lv2_master_printf(" cgm iter 0 rel. residues: ");
    for(int ishift=0;ishift<nshift;ishift++) verbosity_lv2_master_printf("%1.4e  ",1.0);
    verbosity_lv2_master_printf("\n");
    
    int final_iter;
    double final_res[nshift];
    
    int iter=0;
    
    //     -betaa=1
    double betaa=1;
    
    //     -zps=zas=1
    //     -alphas=0
    double zps[nshift],zas[nshift],alphas[nshift];
    double zfs[nshift],betas[nshift];
    int run_flag[nshift],nrun_shift=nshift;
    for(int ishift=0;ishift<nshift;ishift++)
      {
	zps[ishift]=zas[ishift]=1;
	alphas[ishift]=0;
	run_flag[ishift]=1;
      }
    
    //     -alpha=0
    double alpha=0;
    
    //     -rr=(r,r)=source_norm
    double rr=source_norm;
#ifdef CGM_DEBUG
    verbosity_lv3_master_printf("rr: %16.16lg\n",rr);
#endif
    
    double rfrf,pap,betap;
    double res[nshift];
    
    do
      {
	//     this is already iteration 0
	final_iter=(++iter);
	
	//     -s=Ap
	if(use_async_communications && iter>1) CGM_FINISH_COMMUNICATING_BORDERS(p);
	
	if(IS_MASTER_THREAD) cgm_inv_over_time+=take_time();
	APPLY_OPERATOR(s,CGM_OPERATOR_PARAMETERS 0,p);
	if(IS_MASTER_THREAD) cgm_inv_over_time-=take_time();
	
	//     -pap=(p,s)=(p,Ap)
	double_vector_glb_scalar_prod(&pap,(double*)p,(double*)s,BULK_VOL*NDOUBLES_PER_SITE);
#ifdef CGM_DEBUG
	verbosity_lv3_master_printf("pap: %16.16lg (ap[0]: %16.16lg)\n",pap,((double*)s)[0]);
	for(int i=0;i<BULK_VOL*NDOUBLES_PER_SITE;i++)
	  verbosity_lv3_master_printf("%d %lg\n",i,((double*)s)[i]);
#endif
	//     calculate betaa=rr/pap=(r,r)/(p,Ap)
	betap=betaa;
	betaa=-rr/pap;
#ifdef CGM_DEBUG
	verbosity_lv3_master_printf("betap: %16.16lg, betaa: %16.16lg\n",betap,betaa);
#endif
	
	//     calculate
	//     -zfs
	//     -betas
	//     -x
	for(int ishift=0;ishift<nshift;ishift++)
	  {
	    if(run_flag[ishift]==1)
	      {
		double ratio=betap/(betaa*alpha*(1-zas[ishift]/zps[ishift])+betap*(1-shift[ishift]*betaa));
		if(std::isnan(ratio)) crash("nanned");
		zfs[ishift]=zas[ishift]*ratio;
		betas[ishift]=betaa*ratio;
		
#ifdef CGM_DEBUG
		verbosity_lv3_master_printf("ishift %d [%lg] zas: %16.16lg, zps: %16.16lg, "
					    "zfs: %16.16lg, betas: %16.16lg\n",
					    ishift,shift[ishift],zas[ishift],zps[ishift],zfs[ishift],betas[ishift]);
#endif
		double_vector_summ_double_vector_prod_double((double*)(sol[ishift]),(double*)(sol[ishift]),(double*)(ps[ishift]),-betas[ishift],BULK_VOL*NDOUBLES_PER_SITE,DO_NOT_SET_FLAGS);
	      }
	    THREAD_BARRIER();
	  }
	
	//     calculate
	//     -r'=r+betaa*s=r+beta*Ap
	//     -rfrf=(r',r')
	double_vector_summ_double_vector_prod_double((double*)r,(double*)r,(double*)s,betaa,BULK_VOL*NDOUBLES_PER_SITE);
	double_vector_glb_scalar_prod(&rfrf,(double*)r,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
#ifdef CGM_DEBUG
	verbosity_lv3_master_printf("rfrf: %16.16lg\n",rfrf);
#endif
	
	//     calculate alpha=rfrf/rr=(r',r')/(r,r)
	alpha=rfrf/rr;
	
	//     calculate p'=r'+p*alpha
	double_vector_summ_double_vector_prod_double((double*)p,(double*)r,(double*)p,alpha,BULK_VOL*NDOUBLES_PER_SITE);
	
	//start the communications of the border
	if(use_async_communications)  CGM_START_COMMUNICATING_BORDERS(p);
	
	//     calculate 
	//     -alphas=alpha*zfs*betas/zas*beta
	//     -ps'=r'+alpha*ps
	for(int ishift=0;ishift<nshift;ishift++)
	  {
	    if(run_flag[ishift]==1)
	      {
		alphas[ishift]=alpha*zfs[ishift]*betas[ishift]/(zas[ishift]*betaa);
#ifdef CGM_DEBUG
		verbosity_lv3_master_printf("ishift %d alpha: %16.16lg\n",ishift,alphas[ishift]);
#endif
		double_vector_linear_comb((double*)(ps[ishift]),(double*)r,zfs[ishift],(double*)(ps[ishift]),alphas[ishift],BULK_VOL*NDOUBLES_PER_SITE,DO_NOT_SET_FLAGS);
		
		// shift z
		zps[ishift]=zas[ishift];
		zas[ishift]=zfs[ishift];
	      }
	    THREAD_BARRIER();
	  }
	
	//shift rr
	rr=rfrf;
	
	//check over residual
	if(iter%each==0) verbosity_lv2_master_printf(" cgm iter %d rel. residues: ",iter);
	
	for(int ishift=0;ishift<nshift;ishift++)
	  if(run_flag[ishift])
	    {
	      final_res[ishift]=res[ishift]=rr*zfs[ishift]*zfs[ishift]/source_norm;
	      if(iter%each==0)
		{
		  verbosity_lv2_master_printf("%1.4e  ",res[ishift]);
		  verbosity_lv3_master_printf("%16.16lg  ",res[ishift]);
		}
	      
	      if(res[ishift]<inn_req_res[ishift])
		{
		  run_flag[ishift]=0;
		  nrun_shift--;
		}
	    }
	  else
	    if(iter%each==0)
	      verbosity_lv2_master_printf(" * ");
	
	if(iter%each==0)
	  verbosity_lv2_master_printf("\n");
      }
    while(nrun_shift>0 && iter<niter_max);
    
    if(use_async_communications) CGM_FINISH_COMMUNICATING_BORDERS(p);
    
    //print the final true residue
    for(int ishift=0;ishift<nshift;ishift++)
      {
	double res,w_res,weight,max_res;
	APPLY_OPERATOR(s,CGM_OPERATOR_PARAMETERS IN_SHIFT[ishift],sol[ishift]);
	
	double loc_res=0,locw_res=0,locmax_res=0,loc_weight=0;
	NISSA_PARALLEL_LOOP(i,0,BULK_VOL*NDOUBLES_PER_SITE)
	  {
	    double pdiff=((double*)s)[i]-((double*)source)[i];
	    double psol=((double*)sol[ishift])[i];
	    double plain_res=pdiff*pdiff;
	    double point_weight=1/(psol*psol);
	    
	    loc_res+=plain_res;
	    
	    locw_res+=plain_res*point_weight;
	    loc_weight+=point_weight;
	    if(plain_res>locmax_res) locmax_res=plain_res;
	  }
	NISSA_PARALLEL_LOOP_END;
	
	res=glb_reduce_double(loc_res);
	w_res=glb_reduce_double(locw_res);
	weight=glb_reduce_double(loc_weight);
	max_res=glb_max_double(locmax_res);
	
	w_res=w_res/weight;
	
	verbosity_lv2_master_printf(" ishift %d, rel residue true=%lg approx=%lg commanded=%lg weighted=%lg max=%lg\n",
				    ishift,res/source_norm,final_res[ishift],inn_req_res[ishift],w_res,max_res);
	if(res/source_norm>=2*final_res[ishift])
	  master_printf("WARNING: true residue for shift %d (%lg) much larger than expected one (%lg)\n",
			ishift,res/source_norm,final_res[ishift]);
      }
    
    verbosity_lv1_master_printf(" Total cgm iterations: %d\n",final_iter);
    
    //check if not converged
    if(final_iter==niter_max) crash("exit without converging");
    
    for(int ishift=0;ishift<nshift;ishift++) nissa_free(ps[ishift]);
    nissa_free(s);
    nissa_free(p);
    nissa_free(r);
    CGM_ADDITIONAL_VECTORS_FREE();
    
    if(IS_MASTER_THREAD) cgm_inv_over_time+=take_time();
    
#ifdef CG_128_INVERT
    //if 128 bit precision required refine the solution
    if(use_128_bit_precision)
      {
	verbosity_lv1_master_printf("\nRefining the solution in quaduple precision using cg solver\n");
	for(int ishift=0;ishift<nshift;ishift++)
	  CG_128_INVERT(sol[ishift],sol[ishift],CG_128_ADDITIONAL_PARAMETERS_CALL shift[ishift],niter_max,ext_req_res[ishift],source);
      }
#endif
  }
  THREADABLE_FUNCTION_END

//run higher shifts up to common precision
#if CGM_NARG == 0
THREADABLE_FUNCTION_6ARG(CGM_INVERT_RUN_HM_UP_TO_COMM_PREC, BASETYPE**,sol, double*,shift, int,nshift, int,niter_max, double,req_res, BASETYPE*,source)
#elif CGM_NARG == 1
THREADABLE_FUNCTION_7ARG(CGM_INVERT_RUN_HM_UP_TO_COMM_PREC, BASETYPE**,sol, AT1,A1, double*,shift, int,nshift, int,niter_max, double,req_res, BASETYPE*,source)
#elif CGM_NARG == 2
THREADABLE_FUNCTION_8ARG(CGM_INVERT_RUN_HM_UP_TO_COMM_PREC, BASETYPE**,sol, AT1,A1, AT2,A2, double*,shift, int,nshift, int,niter_max, double,req_res, BASETYPE*,source)
#elif CGM_NARG == 3
THREADABLE_FUNCTION_9ARG(CGM_INVERT_RUN_HM_UP_TO_COMM_PREC, BASETYPE**,sol, AT1,A1, AT2,A2, AT3,A3, double*,shift, int,nshift, int,niter_max, double,req_res, BASETYPE*,source)
#elif CGM_NARG == 4
THREADABLE_FUNCTION_10ARG(CGM_INVERT_RUN_HM_UP_TO_COMM_PREC, BASETYPE**,sol, AT1,A1, AT2,A2, AT3,A3, AT4,A4, double*,shift, int,nshift, int,niter_max, double,req_res, BASETYPE*,source)
#elif CGM_NARG == 5
THREADABLE_FUNCTION_11ARG(CGM_INVERT_RUN_HM_UP_TO_COMM_PREC, BASETYPE**,sol, AT1,A1, AT2,A2, AT3,A3, AT4,A4, AT5,A5, double*,shift, int,nshift, int,niter_max, double,req_res, BASETYPE*,source)
#endif
  {
    double req_res_int[nshift];
    for(int ishift=0;ishift<nshift;ishift++) req_res_int[ishift]=req_res;
    CGM_INVERT(sol,CGM_ADDITIONAL_PARAMETERS_CALL shift,nshift,niter_max,req_res_int,source);
  }
  THREADABLE_FUNCTION_END
  
  //put outside to be common with single version
#include "cgm_invert_template_summsol_threaded.hpp"
}

#undef BASETYPE
#undef NDOUBLES_PER_SITE
#undef BULK_VOL
#undef BORD_VOL

#undef APPLY_OPERATOR

#undef CGM_INVERT_RUN_HM_TO_COMM_PREC
#undef SUMM_SRC_AND_ALL_INV_CGM
#undef CGM_INVERT
#undef CGM_START_COMMUNICATING_BORDERS
#undef CGM_FINISH_COMMUNICATING_BORDERS
#undef CGM_NPOSSIBLE_REQUEST

#undef CGM_ADDITIONAL_VECTORS_ALLOCATION
#undef CGM_ADDITIONAL_VECTORS_FREE

#undef CGM_ADDITIONAL_PARAMETERS_CALL
#undef CGM_ADDITIONAL_PARAMETERS_PROTO
#undef CGM_OPERATOR_PARAMETERS

#undef CG_128_INVERT

#undef CGM_NARG
