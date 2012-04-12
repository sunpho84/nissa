#pragma once

void cgmm2s_invert(basetype **sol,quad_su3 **conf,double *m2,int nmass,int niter,double st_res,double st_minres,int st_crit,basetype *source)
{
  const int each=10;
  
  basetype *t=nissa_malloc("DD_temp",bulk_vol+bord_vol,basetype);
  basetype *s=nissa_malloc("s",bulk_vol,basetype);
  basetype *r=nissa_malloc("r",bulk_vol,basetype);
  basetype *p=nissa_malloc("p",bulk_vol+bord_vol,basetype);
  basetype *ps[nmass];
  for(int imass=0;imass<nmass;imass++) ps[imass]=nissa_malloc("ps",bulk_vol,basetype);
  
  //sol[*]=0
  int run_flag[nmass],nrun_mass=nmass;
  for(int imass=0;imass<nmass;imass++)
    {
      double_vector_init_to_zero((double*)(sol[imass]),bulk_vol*ndoubles_per_site);
      run_flag[imass]=1;
    }
  
  //     -p=source
  //     -r=source
  //     -calculate Rr=(r,r)
  double_vector_copy((double*)p,(double*)source,bulk_vol*ndoubles_per_site);
  double_vector_copy((double*)r,(double*)source,bulk_vol*ndoubles_per_site);
  double rr=double_vector_glb_scalar_prod((double*)r,(double*)r,bulk_vol*ndoubles_per_site);
  
  //writes source norm
  double source_norm=rr;
  verbosity_lv2_master_printf(" Source norm: %lg\n",source_norm);
  if(source_norm==0 || isnan(source_norm)) crash("invalid norm: %lg",source_norm);
  
  //writes initial resiude
  verbosity_lv2_master_printf(" cgmms iter 0 rel. residues: ");
  for(int imass=0;imass<nmass;imass++) verbosity_lv2_master_printf("%1.4e  ",1.0);
  verbosity_lv2_master_printf("\n");
  
  //     -betaa=1
  double betaa=1;
  
  //     -ps=source
  //     -zps=zas=1
  //     -alphas=0
  double zps[nmass],zas[nmass],alphas[nmass];
  for(int imass=0;imass<nmass;imass++)
    {
      double_vector_copy((double*)(ps[imass]),(double*)source,bulk_vol*ndoubles_per_site);
      
      zps[imass]=zas[imass]=1;
      alphas[imass]=0;
    }
  
  //     -alpha=0
  double alpha=0;
      
  int iter=0;
  double final_res[nmass];
  int nrequest=0;
  MPI_Request request[cgmm2s_npossible_requests];
  do
    {
      double zfs[nmass],betas[nmass];
      
      //     -s=Ap
      //if(nrequest!=0) finish_communicating_ev_color_borders(&nrequest,request,p);
      apply_offdiagonal_operator(s,conf,t,p);
      
      //     -pap=(p,s)=(p,Ap)
      double pap=double_vector_glb_scalar_prod((double*)p,(double*)s,bulk_vol*ndoubles_per_site);
      
      //     calculate betaa=rr/pap=(r,r)/(p,Ap)
      double betap=betaa;
      betaa=-rr/pap;
      
      //     calculate 
      //     -zfs
      //     -betas
      //     -x
      for(int imass=0;imass<nmass;imass++)
	{
	  if(run_flag[imass]==1)
	    {
	      zfs[imass]=zas[imass]*betap/(betaa*alpha*(1-zas[imass]/zps[imass])+betap*(1-m2[imass]*betaa));
	      betas[imass]=betaa*zfs[imass]/zas[imass];
	      
	      double_vector_subt_double_vector_prod_double((double*)(sol[imass]),(double*)(sol[imass]),(double*)(ps[imass]),betas[imass],bulk_vol*ndoubles_per_site);
	    }
	}

      //     calculate
      //     -r'=r+betaa*s=r+beta*Ap
      //     -rfrf=(r',r')
      double_vector_summ_double_vector_prod_double((double*)r,(double*)r,(double*)s,betaa,bulk_vol*ndoubles_per_site);
      double rfrf=double_vector_glb_scalar_prod((double*)r,(double*)r,bulk_vol*ndoubles_per_site);
      
      //     calculate alpha=rfrf/rr=(r',r')/(r,r)
      alpha=rfrf/rr;
      
      //     calculate p'=r'+p*alpha
      double_vector_summ_double_vector_prod_double((double*)p,(double*)r,(double*)p,alpha,bulk_vol*ndoubles_per_site);
      
      //start the communications of the border
      //cgmm2s_start_communicating_borders(&nrequest,request,p);
      
      //     calculate 
      //     -alphas=alpha*zfs*betas/zas*beta
      //     -ps'=r'+alpha*ps
      for(int imass=0;imass<nmass;imass++)
	if(run_flag[imass]==1)
	  {
	    alphas[imass]=alpha*zfs[imass]*betas[imass]/(zas[imass]*betaa);
	    double_vector_linear_comb((double*)(ps[imass]),(double*)r,zfs[imass],(double*)(ps[imass]),alphas[imass],bulk_vol*ndoubles_per_site);

	    // shift z
	    zps[imass]=zas[imass];
	    zas[imass]=zfs[imass];
	  }
      
      //shift rr
      rr=rfrf;
      
      iter++;
      
      //check over residual
      if(iter%each==0) verbosity_lv2_master_printf(" cgmms iter %d rel. residues: ",iter);
      for(int imass=0;imass<nmass;imass++)
	if(run_flag[imass])
	  {
	    int fini=0;
	    if(st_crit==sc_standard || st_crit==sc_unilevel)
	      {
		final_res[imass]=rr*zfs[imass]*zfs[imass]/source_norm;
		if(st_crit==sc_standard) fini=(final_res[imass]<st_minres || final_res[0]<st_res);
		else fini=final_res[imass]<st_res;
		if(iter%each==0) verbosity_lv2_master_printf("%1.4e  ",final_res[imass]);
	      }
	    else crash("unkwnown stopping criterion");
	    
	    if(fini)
	      {
		run_flag[imass]=0;
		nrun_mass--;
	      }
	  } 
      if(iter%each==0) verbosity_lv2_master_printf("\n");
    }
  while(nrun_mass>0 && iter<niter);
  
  //print the final true residue
  
  for(int imass=0;imass<nmass;imass++)
    {
      double res,w_res,weight,max_res;
      apply_full_operator(s,conf,t,sqrt(m2[imass]),sol[imass]);
      {
	double loc_res=0;
	double locw_res=0;
	double locmax_res=0;
	double loc_weight=0;
	
	complex *ds=(complex*)s,*dsour=(complex*)source,*dsol=(complex*)(sol[imass]);
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
	
	verbosity_lv1_master_printf(" imass %d, rel residue true=%g approx=%g weighted=%g max=%g\n",imass,res/source_norm,final_res[imass],w_res,max_res);
      }
    }  
  
  verbosity_lv1_master_printf(" Total cgmms iterations: %d\n",iter);
  
  for(int imass=0;imass<nmass;imass++) nissa_free(ps[imass]);
  nissa_free(s);
  nissa_free(p);
  nissa_free(r);
  nissa_free(t);
}

//return all the masses summed together
void summ_src_and_all_inv_cgmm2s(basetype *sol,quad_su3 **conf,rat_approx *appr,int niter,double st_res,double st_minres,int st_crit,basetype *source)
{
  //allocate temporary single solutions
  basetype *temp[appr->nterms];
  for(int iterm=0;iterm<appr->nterms;iterm++)
    temp[iterm]=nissa_malloc("temp",bulk_vol+bord_vol,basetype);
  
  //call multi-mass solver
  cgmm2s_invert(temp,conf,appr->poles,appr->nterms,niter,st_res,st_minres,st_crit,source);
  
  //summ all the masses
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

#undef apply_offdiagonal_operator
#undef apply_full_operator

#undef summ_src_and_all_inv_cgmm2s
#undef cgmm2s_invert
#undef cgmm2s_start_communicating_borders
#undef cgmm2s_npossible_requests
