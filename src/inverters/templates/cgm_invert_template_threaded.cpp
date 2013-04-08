#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <omp.h>
 #include <unistd.h>

#include "../../routines/ios.h"
#include "../../routines/openmp.h"
#include "../../routines/mpi.h"
#include "../../base/openmp_macros.h"

extern double cgm_inv_over_time;
extern int ncgm_inv;

/*
  This is the prorotipe for a multi-shift inverter.
  The calls to the operator, the internal vectors definitions and the additional parameters must be defined thorugh macro.
  The file must be included inside another file defining all the macros.
  See "cgm_invert_tmQ2.c" as an example.
*/

#if CGM_NARG >= 4
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
#endif
{
  GET_THREAD_ID();
  
#ifdef CG_128_INVERT
  //limit inner solver precision
  double max_inner_solver=1.0e-25;
  //used for inner solver in the case of 128 bit precision
  double inn_req_res[nshift];
  for(int ishift=0;ishift<nshift;ishift++)
    if(nissa_use_128_bit_precision && ext_req_res[ishift]<max_inner_solver)
      {
	verbosity_lv2_master_printf("changing the inner solver residue for shift %d to %lg\n",ishift,max_inner_solver);
	inn_req_res[ishift]=max_inner_solver;
      }
    else inn_req_res[ishift]=ext_req_res[ishift];
#else
  double *inn_req_res=ext_req_res;
#endif
	
  if(IS_MASTER_THREAD)
    {
      ncgm_inv++;
      cgm_inv_over_time-=take_time();
    }
  
  int each_list[4]={0,100,10,1},each;
  if(nissa_verbosity>=3) each=1;
  else each=each_list[nissa_verbosity];
  
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
  if(source_norm==0 || isnan(source_norm)) crash("invalid norm: %lg",source_norm);
  
  //writes initial residue
  verbosity_lv2_master_printf(" cgm iter 0 rel. residues: ");
  for(int ishift=0;ishift<nshift;ishift++) verbosity_lv2_master_printf("%1.4e  ",1.0);
  verbosity_lv2_master_printf("\n");
  
  int nrequest=0;
  int final_iter;
  MPI_Request request[CGM_NPOSSIBLE_REQUESTS];
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
  
  double rfrf,pap,betap;
  double res[nshift];
  
  do
    {
      //     this is already iteration 0
      final_iter=(++iter);
      
      //     -s=Ap
      if(nissa_use_async_communications) CGM_FINISH_COMMUNICATING_BORDERS(&nrequest,request,p);
      if(IS_MASTER_THREAD) cgm_inv_over_time+=take_time();
      APPLY_OPERATOR(s,CGM_OPERATOR_PARAMETERS shift[0],p);
      if(IS_MASTER_THREAD) cgm_inv_over_time-=take_time();
      
      //     -pap=(p,s)=(p,Ap)
      double_vector_glb_scalar_prod(&pap,(double*)p,(double*)s,BULK_VOL*NDOUBLES_PER_SITE);
      
      //     calculate betaa=rr/pap=(r,r)/(p,Ap)
      betap=betaa;
      betaa=-rr/pap;
      
      //     calculate 
      //     -zfs
      //     -betas
      //     -x
      for(int ishift=0;ishift<nshift;ishift++)
	{
	  if(run_flag[ishift]==1)
	    {
	      zfs[ishift]=zas[ishift]*betap/(betaa*alpha*(1-zas[ishift]/zps[ishift])+betap*(1-(shift[ishift]-shift[0])*betaa));
	      betas[ishift]=betaa*zfs[ishift]/zas[ishift];
	      
	      double_vector_summ_double_vector_prod_double((double*)(sol[ishift]),(double*)(sol[ishift]),(double*)(ps[ishift]),-betas[ishift],BULK_VOL*NDOUBLES_PER_SITE,DO_NOT_SET_FLAGS);
	    }
	}
	
      //     calculate
      //     -r'=r+betaa*s=r+beta*Ap
      //     -rfrf=(r',r')
      double_vector_summ_double_vector_prod_double((double*)r,(double*)r,(double*)s,betaa,BULK_VOL*NDOUBLES_PER_SITE,DO_NOT_SET_FLAGS);
      double_vector_glb_scalar_prod(&rfrf,(double*)r,(double*)r,BULK_VOL*NDOUBLES_PER_SITE);
      
      //     calculate alpha=rfrf/rr=(r',r')/(r,r)
      alpha=rfrf/rr;
	
      //     calculate p'=r'+p*alpha
      double_vector_summ_double_vector_prod_double((double*)p,(double*)r,(double*)p,alpha,BULK_VOL*NDOUBLES_PER_SITE);
      
      //start the communications of the border
      if(nissa_use_async_communications) CGM_START_COMMUNICATING_BORDERS(&nrequest,request,p);
      
      //     calculate 
      //     -alphas=alpha*zfs*betas/zas*beta
      //     -ps'=r'+alpha*ps
      for(int ishift=0;ishift<nshift;ishift++)
	if(run_flag[ishift]==1)
	  {
	    alphas[ishift]=alpha*zfs[ishift]*betas[ishift]/(zas[ishift]*betaa);

	    double_vector_linear_comb((double*)(ps[ishift]),(double*)r,zfs[ishift],(double*)(ps[ishift]),alphas[ishift],BULK_VOL*NDOUBLES_PER_SITE,DO_NOT_SET_FLAGS);
	     
	    // shift z
	    zps[ishift]=zas[ishift];
	    zas[ishift]=zfs[ishift];
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
		if(nissa_verbosity==2) master_printf("%1.4e  ",res[ishift]);
		if(nissa_verbosity==3) master_printf("%16.16lg  ",res[ishift]);
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
  
  if(nissa_use_async_communications) CGM_FINISH_COMMUNICATING_BORDERS(&nrequest,request,p);
  
  //print the final true residue
  for(int ishift=0;ishift<nshift;ishift++)
    {
      double res,w_res,weight,max_res;
      APPLY_OPERATOR(s,CGM_OPERATOR_PARAMETERS shift[ishift],sol[ishift]);

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
      
      res=glb_reduce_double(loc_res);
      w_res=glb_reduce_double(locw_res);
      weight=glb_reduce_double(loc_weight);
      //max_res=glb_max_double(locmax_res); //to be fixed

      w_res=w_res/weight;
      
      verbosity_lv2_master_printf(" ishift %d, rel residue true=%lg approx=%lg commanded=%lg weighted=%lg max=%lg\n",
				  ishift,res/source_norm,final_res[ishift],inn_req_res[ishift],w_res,locmax_res);
    }  
  
  verbosity_lv1_master_printf(" Total cgm iterations: %d\n",final_iter);
  
  for(int ishift=0;ishift<nshift;ishift++) nissa_free(ps[ishift]);
  nissa_free(s);
  nissa_free(p);
  nissa_free(r);
  CGM_ADDITIONAL_VECTORS_FREE();
  
  if(IS_MASTER_THREAD) cgm_inv_over_time+=take_time();
  
#ifdef CG_128_INVERT
  //if 128 bit precision required refine the solution
  if(nissa_use_128_bit_precision)
    {
      verbosity_lv1_master_printf("\nRefining the solution in quaduple precision using cg solver\n");
      for(int ishift=0;ishift<nshift;ishift++)
	CG_128_INVERT(sol[ishift],sol[ishift],CG_128_ADDITIONAL_PARAMETERS_CALL shift[ishift],niter_max,5,ext_req_res[ishift],source);
    }
#endif
}}

//run higher shifts up to common precision
#if CGM_NARG == 0
THREADABLE_FUNCTION_6ARG(CGM_INVERT_RUN_HM_UP_TO_COMM_PREC, BASETYPE**,sol, double*,shift, int,nshift, int,niter_max, double,req_res, BASETYPE*,source)
#elif CGM_NARG == 1
THREADABLE_FUNCTION_7ARG(CGM_INVERT_RUN_HM_UP_TO_COMM_PREC, BASETYPE**,sol, AT1,A1, double*,shift, int,nshift, int,niter_max, double,req_res, BASETYPE*,source)
#elif CGM_NARG == 2
THREADABLE_FUNCTION_8ARG(CGM_INVERT_RUN_HM_UP_TO_COMM_PREC, BASETYPE**,sol, AT1,A1, AT2,A2, double*,shift, int,nshift, int,niter_max, double,req_res, BASETYPE*,source)
#elif CGM_NARG == 3
THREADABLE_FUNCTION_9ARG(CGM_INVERT_RUN_HM_UP_TO_COMM_PREC, BASETYPE**,sol, AT1,A1, AT2,A2, AT3,A3, double*,shift, int,nshift, int,niter_max, double,req_res, BASETYPE*,source)
#endif
{
  double req_res_int[nshift];
  for(int ishift=0;ishift<nshift;ishift++) req_res_int[ishift]=req_res;
  CGM_INVERT(sol,CGM_ADDITIONAL_PARAMETERS_CALL shift,nshift,niter_max,req_res_int,source);
}}

//return all the shifts summed together
#if CGM_NARG == 0
THREADABLE_FUNCTION_5ARG(SUMM_SRC_AND_ALL_INV_CGM, BASETYPE*,sol, rat_approx_t*,appr, int,niter_max, double,req_res, BASETYPE*,source)
#elif CGM_NARG == 1
THREADABLE_FUNCTION_6ARG(SUMM_SRC_AND_ALL_INV_CGM, BASETYPE*,sol, AT1,A1, rat_approx_t*,appr, int,niter_max, double,req_res, BASETYPE*,source)
#elif CGM_NARG == 2
THREADABLE_FUNCTION_7ARG(SUMM_SRC_AND_ALL_INV_CGM, BASETYPE*,sol, AT1,A1, AT2,A2, rat_approx_t*,appr, int,niter_max, double,req_res, BASETYPE*,source)
#elif CGM_NARG == 3
THREADABLE_FUNCTION_8ARG(SUMM_SRC_AND_ALL_INV_CGM, BASETYPE*,sol, AT1,A1, AT2,A2, AT3,A3, rat_approx_t*,appr, int,niter_max, double,req_res, BASETYPE*,source)
#endif
{
  GET_THREAD_ID();
  
  //allocate temporary single solutions
  BASETYPE *temp[appr->degree];
  for(int iterm=0;iterm<appr->degree;iterm++)
    temp[iterm]=nissa_malloc(combine("temp%d",iterm).c_str(),BULK_VOL+BORD_VOL,BASETYPE);
  
  //call multi-shift solver
  CGM_INVERT_RUN_HM_UP_TO_COMM_PREC(temp,CGM_ADDITIONAL_PARAMETERS_CALL appr->poles,appr->degree,niter_max,req_res,source);
  
  //summ all the shifts
  NISSA_PARALLEL_LOOP(i,0,BULK_VOL*NDOUBLES_PER_SITE)
    {
      ((double*)sol)[i]=appr->cons*((double*)source)[i];
      for(int iterm=0;iterm<appr->degree;iterm++)
	((double*)sol)[i]+=appr->weights[iterm]*((double*)(temp[iterm]))[i];
    }
  
  set_borders_invalid(sol);
  
  //free temp vectors
  for(int iterm=0;iterm<appr->degree;iterm++)
    nissa_free(temp[iterm]);
}}

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
