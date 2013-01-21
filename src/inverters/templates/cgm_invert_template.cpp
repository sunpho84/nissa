/*
  This is the prorotipe for a multi-shift inverter.
  The calls to the operator, the internal vectors definitions and the additional parameters must be defined thorugh macro.
  The file must be included inside another file defining all the macros.
  See "cgm_invert_tmQ2.c" as an example.
*/

void cgm_invert(basetype **sol,cgm_additional_parameters_proto,double *shift,int nshift,int niter_max,double *req_res,basetype *source)
{
  const int each=10;
  
  //macro to be defined externally, allocating all the required additional vectors
  cgm_additional_vectors_allocation();
  
  basetype *s=nissa_malloc("s",bulk_vol,basetype);
  basetype *r=nissa_malloc("r",bulk_vol,basetype);
  basetype *p=nissa_malloc("p",bulk_vol+bord_vol,basetype);
  basetype *ps[nshift];
  for(int ishift=0;ishift<nshift;ishift++) ps[ishift]=nissa_malloc("ps",bulk_vol,basetype);
  
  //     -sol[*]=0
  //     -ps[*]=source
  for(int ishift=0;ishift<nshift;ishift++)
    {
      double_vector_copy((double*)(ps[ishift]),(double*)source,bulk_vol*ndoubles_per_site);
      double_vector_init_to_zero((double*)(sol[ishift]),bulk_vol*ndoubles_per_site);
    }
  
  //     -p=source
  //     -r=source
  //     -calculate source_norm=(r,r)
  double_vector_copy((double*)p,(double*)source,bulk_vol*ndoubles_per_site);
  double_vector_copy((double*)r,(double*)source,bulk_vol*ndoubles_per_site);
  double source_norm=double_vector_glb_scalar_prod((double*)r,(double*)r,bulk_vol*ndoubles_per_site);
  
  //writes source norm
  verbosity_lv2_master_printf(" Source norm: %lg\n",source_norm);
  if(source_norm==0 || isnan(source_norm)) crash("invalid norm: %lg",source_norm);
  
  //writes initial resiude
  verbosity_lv2_master_printf(" cgm iter 0 rel. residues: ");
  for(int ishift=0;ishift<nshift;ishift++) verbosity_lv2_master_printf("%1.4e  ",1.0);
  verbosity_lv2_master_printf("\n");
  
  int nrequest=0;
  int iter=0;
  MPI_Request request[cgm_npossible_requests];
  double final_res[nshift];
  
#pragma omp parallel
  {
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
#pragma omp single
	iter++;
	
	//     -s=Ap
#pragma omp single
	if(nissa_use_async_communications && nrequest!=0) cgm_finish_communicating_borders(nrequest,request,p);
	
	apply_operator(s,cgm_operator_parameters,shift[0],p);
	
	//     -pap=(p,s)=(p,Ap)
	pap=double_vector_glb_scalar_prod((double*)p,(double*)s,bulk_vol*ndoubles_per_site);
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
		
		double_vector_subt_double_vector_prod_double((double*)(sol[ishift]),(double*)(sol[ishift]),(double*)(ps[ishift]),betas[ishift],bulk_vol*ndoubles_per_site);
	      }
	  }
	
	//     calculate
	//     -r'=r+betaa*s=r+beta*Ap
	//     -rfrf=(r',r')
	double_vector_summ_double_vector_prod_double((double*)r,(double*)r,(double*)s,betaa,bulk_vol*ndoubles_per_site);
	rfrf=double_vector_glb_scalar_prod((double*)r,(double*)r,bulk_vol*ndoubles_per_site);
	
	//     calculate alpha=rfrf/rr=(r',r')/(r,r)
	alpha=rfrf/rr;
	
	//     calculate p'=r'+p*alpha
	double_vector_summ_double_vector_prod_double((double*)p,(double*)r,(double*)p,alpha,bulk_vol*ndoubles_per_site);
	
	//start the communications of the border
#pragma omp single
	if(nissa_use_async_communications) cgm_start_communicating_borders(nrequest,request,p);
	
	//     calculate 
	//     -alphas=alpha*zfs*betas/zas*beta
	//     -ps'=r'+alpha*ps
	for(int ishift=0;ishift<nshift;ishift++)
	  if(run_flag[ishift]==1)
	   {
	    alphas[ishift]=alpha*zfs[ishift]*betas[ishift]/(zas[ishift]*betaa);
	    
	    double_vector_linear_comb((double*)(ps[ishift]),(double*)r,zfs[ishift],(double*)(ps[ishift]),alphas[ishift],bulk_vol*ndoubles_per_site);
	     
	    // shift z
	    zps[ishift]=zas[ishift];
	    zas[ishift]=zfs[ishift];
	   }
	
	//shift rr
	rr=rfrf;
	
	//check over residual
#pragma omp single
	if(iter%each==0) verbosity_lv2_master_printf(" cgm iter %d rel. residues: ",iter);
	for(int ishift=0;ishift<nshift;ishift++)
	  if(run_flag[ishift])
	    {
	      final_res[ishift]=res[ishift]=rr*zfs[ishift]*zfs[ishift]/source_norm;
#pragma omp single
	      if(iter%each==0) verbosity_lv2_master_printf("%1.4e  ",res[ishift]);
	      
	      if(res[ishift]<req_res[ishift])
		{
		  run_flag[ishift]=0;
		  nrun_shift--;
		}
	    }
	  else
#pragma omp single
	    if(iter%each==0) verbosity_lv2_master_printf(" * ");
#pragma omp single
	if(iter%each==0) verbosity_lv2_master_printf("\n");
      }
    while(nrun_shift>0 && iter<niter_max);
  }
  
  if(nissa_use_async_communications && nrequest!=0) cgm_finish_communicating_borders(nrequest,request,p);
  
  //print the final true residue
  for(int ishift=0;ishift<nshift;ishift++)
    {
      double res,w_res,weight,max_res;
      apply_operator(s,cgm_operator_parameters,shift[ishift],sol[ishift]);
      {
	double loc_res=0;
	double locw_res=0;
	double locmax_res=0;
	double loc_weight=0;
	
	complex *ds=(complex*)s,*dsour=(complex*)source,*dsol=(complex*)(sol[ishift]);
	for(int i=0;i<bulk_vol*ndoubles_per_site/2;i++)
	  {
	    ds[i][0]-=dsour[i][0];
	    ds[i][1]-=dsour[i][1];
	    double plain_res=ds[i][0]*ds[i][0]+ds[i][1]*ds[i][1];
	    double point_weight=1/(dsol[i][0]*dsol[i][0]+dsol[i][1]*dsol[i][1]);
	    
	    loc_res+=plain_res;
	    
	    locw_res+=plain_res*point_weight;
	    loc_weight+=point_weight;
	    if(plain_res>locmax_res) locmax_res=plain_res;
	  }
	set_borders_invalid(s);
	
	MPI_Reduce(&loc_res,&res,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&locw_res,&w_res,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&loc_weight,&weight,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&locmax_res,&max_res,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	
	w_res=w_res/weight;
	
	verbosity_lv2_master_printf(" ishift %d, rel residue true=%lg approx=%lg commanded=%lg weighted=%lg max=%lg\n",
				    ishift,res/source_norm,final_res[ishift],req_res[ishift],w_res,max_res);
      }
    }  
  
  verbosity_lv1_master_printf(" Total cgm iterations: %d\n",iter);
  
  for(int ishift=0;ishift<nshift;ishift++) nissa_free(ps[ishift]);
  nissa_free(s);
  nissa_free(p);
  nissa_free(r);
  cgm_additional_vectors_free();
  
#ifdef cg_128_invert
  //if 128 bit precision required refine the solution
  if(nissa_use_128_bit_precision)
    {
      verbosity_lv1_master_printf("\nRefining the solution in quaduple precision using cg solver\n");
      for(int ishift=0;ishift<nshift;ishift++)
	cg_128_invert(sol[ishift],sol[ishift],cg_128_additional_parameters_call,shift[ishift],niter_max,5,req_res[ishift],source);
    }
#endif
}

//run higher shifts up to common precision
void cgm_invert_run_hm_up_to_comm_prec(basetype **sol,cgm_additional_parameters_proto,double *shift,int nshift,int niter_max,double req_res,basetype *source)
{
  double req_res_int[nshift];
  for(int ishift=0;ishift<nshift;ishift++) req_res_int[ishift]=req_res;
  cgm_invert(sol,cgm_additional_parameters_call,shift,nshift,niter_max,req_res_int,source);
}

//return all the shifts summed together
void summ_src_and_all_inv_cgm(basetype *sol,cgm_additional_parameters_proto,rat_approx_type *appr,int niter_max,double req_res,basetype *source)
{
  //allocate temporary single solutions
  basetype *temp[appr->nterms];
  for(int iterm=0;iterm<appr->nterms;iterm++)
    temp[iterm]=nissa_malloc("temp",bulk_vol+bord_vol,basetype);
  
  //call multi-shift solver
  cgm_invert_run_hm_up_to_comm_prec(temp,cgm_additional_parameters_call,appr->poles,appr->nterms,niter_max,req_res,source);
  
  //summ all the shifts
  double *dsol=(double*)sol,*dsource=(double*)source;
  for(int i=0;i<bulk_vol*ndoubles_per_site;i++)
    {
      dsol[i]=appr->cons*dsource[i];
      for(int iterm=0;iterm<appr->nterms;iterm++)
	dsol[i]+=appr->weights[iterm]*(((double*)(temp[iterm]))[i]);
    }
  
  set_borders_invalid(sol);
  
  //free temp vectors
  for(int iterm=0;iterm<appr->nterms;iterm++)
    nissa_free(temp[iterm]);
}

#undef basetype
#undef ndoubles_per_site
#undef bulk_vol
#undef bord_vol

#undef apply_operator

#undef cgm_invert_run_hm_up_to_comm_prec
#undef summ_src_and_all_inv_cgm
#undef cgm_invert
#undef cgm_start_communicating_borders
#undef cgm_finish_communicating_borders
#undef cgm_npossible_requests

#undef cgm_additional_vectors_allocation
#undef cgm_additional_vectors_free

#undef cgm_additional_parameters_call
#undef cgm_additional_parameters_proto
#undef cgm_operator_parameters

#undef cg_128_invert
