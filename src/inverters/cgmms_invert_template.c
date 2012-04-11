#pragma once

double cgmm2s_calculate_weighted_residue(basetype *source,quad_su3 **conf,double m2,basetype *s,basetype *t,int dinf,basetype *sol)
{
  apply_full_operator(s,conf,t,sqrt(m2),sol);
  
  double loc_weighted_residue=0,tot_weighted_residue;
  double loc_weight=0,tot_weight;
  double *ds=(double*)s,*dsource=(double*)source,*dsol=(double*)sol;
  
  for(int i=0;i<bulk_vol*ndoubles_per_site;i++)
    {
      double nr=(*ds)-(*dsource);
      double ni=(*(ds+1))-(*(dsource+1));
      
      double weight=1/((*dsol)*(*dsol)+(*dsol+1)*(*dsol+1));
      double contrib=(nr*nr+ni*ni)*weight;
      
      if(dinf==2)
        {
          loc_weighted_residue+=contrib;
          loc_weight+=weight;
        }
      else if(contrib>loc_weighted_residue) loc_weighted_residue=contrib;
        
      ds+=2;dsource+=2;dsol+=2;
    }
  
    if(dinf==2)
      {
        MPI_Allreduce(&loc_weighted_residue,&tot_weighted_residue,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&loc_weight,&tot_weight,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
    else MPI_Allreduce(&loc_weighted_residue,&tot_weighted_residue,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  
  if(dinf==2) tot_weighted_residue/=bulk_vol*ndoubles_per_site*tot_weight;
  
  return tot_weighted_residue;
}

int cgmm2s_check_residue(int *run_flag,double *residue_mass,int nrun,double rr,double *zfs,int st_crit,double st_res,double st_res2,int iter,int nmass,double *m2,basetype *source,quad_su3 **conf,basetype *s,basetype *t,double source_norm,basetype **sol)
{
  const int each=10;

  for(int imass=0;imass<nmass;imass++)
    if(run_flag[imass])
      {
	int fini=0;

	if(st_crit==sc_standard||st_crit==sc_unilevel)
	  {
	    residue_mass[imass]=rr*zfs[imass]*zfs[imass]/source_norm;
	    
	    if(st_crit==sc_standard) fini=(residue_mass[imass]<st_res2||residue_mass[0]<st_res);
	    else fini=residue_mass[imass]<st_res;
	  }
	else if(st_crit==sc_weighted_norm2||st_crit==sc_weighted_norm_inf)
	  if(iter%each==0)
	    {  //locally weighted norm
	      if(st_crit==sc_weighted_norm2)
		residue_mass[imass]=cgmm2s_calculate_weighted_residue(source,conf,sqrt(m2[imass]),s,t,2,sol[imass]);
	      else residue_mass[imass]=cgmm2s_calculate_weighted_residue(source,conf,sqrt(m2[imass]),s,t,-1,sol[imass]);
	    }
	
	if(fini)
	  {
	    run_flag[imass]=0;
	    nrun--;
	  }
      } 
  
  if((iter%each==0) && (nissa_inverter_verbosity>=2))
    {    
      master_printf(" cgmms iter %d rel. residues: ",iter);
      for(int imass=0;imass<nmass;imass++)
	if(run_flag[imass])
	  master_printf("%1.4e  ",residue_mass[imass]);
      master_printf("\n");
    }
  
  return nrun;
}

//     -p=source
//     -r=source
//     -calculate Rr=(r,r)
double cgmm2s_init_vector_and_compute_init_residue(basetype *p,basetype *r,basetype *source)
{
  double loc_rr=0;
  double *dsource=(double*)source,*dp=(double*)p,*dr=(double*)r;
  for(int i=0;i<bulk_vol*ndoubles_per_site;i++)
    {
      (*dp)=(*dsource);
      (*dr)=(*dsource);
      loc_rr+=(*dr)*(*dr);
      
      dsource++;dp++;dr++;
    }
  set_borders_invalid(p);
  
  double rr;
  MPI_Allreduce(&loc_rr,&rr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return rr;
}

//     -ps=source
//     -zps=zas=1
//     -alphas=0
void cgmm2s_init_additional_vectors_and_weights(basetype **ps,basetype *source,double *zps,double *zas,double *alphas,int nmass)
{
  for(int imass=0;imass<nmass;imass++)
    {
      double *dps=(double*)(ps[imass]),*dsource=(double*)source;
      for(int i=0;i<bulk_vol*ndoubles_per_site;i++)
	{
	  (*dps)=(*dsource);
	  dps++;dsource++;
	}
      
      zps[imass]=zas[imass]=1;
      alphas[imass]=0;
    }
}

//     -pap=(p,s)=(p,Ap)
double cgmm2s_compute_scalar_product(basetype *p,basetype *s)
{
  double loc_pap=0;
  double *dp=(double*)p,*ds=(double*)s;
  int n=0;
  for(int i=0;i<bulk_vol*ndoubles_per_site;i++)
    {
      loc_pap+=(*dp)*(*ds);
      dp++;ds++;
      n++;
    }

  double pap;
  MPI_Allreduce(&loc_pap,&pap,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return pap;
}

//     -r'=r+betaa*s=r+beta*Ap
//     -rfrf=(r',r')
double cgmm2s_compute_new_residue(basetype *r,basetype *s,double betaa)
{
  double loc_rfrf=0;
  
  double *dr=(double*)r,*ds=(double*)s;
  for(int i=0;i<bulk_vol*ndoubles_per_site;i++)
    {
      (*dr)+=(*ds)*betaa;
      loc_rfrf+=(*dr)*(*dr);
      dr++;ds++;
    }
  
  double rfrf;
  MPI_Allreduce(&loc_rfrf,&rfrf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return rfrf;
}

//     -r'=r+betaa*s=r+beta*Ap
//     -rfrf=(r',r')
void cgmm2s_compute_new_vector(basetype *p,basetype *r,double alpha)
{
  double *dp=(double*)p,*dr=(double*)r;		
  for(int i=0;i<bulk_vol*ndoubles_per_site;i++)
    {
      (*dp)=(*dr)+alpha*(*dp);
      dp++;dr++;
    }
  set_borders_invalid(p);
}

//     -alphas=alpha*zfs*betas/zas*beta
//     -ps'=r'+alpha*ps
void cgmm2s_scale_additional_vectors(basetype **ps,basetype *r,double alpha,double *alphas,double betaa,double *betas,double *zfs,double *zas,int *run_flag,int nmass)
{
  for(int imass=0;imass<nmass;imass++)
    if(run_flag[imass]==1)
      {
	alphas[imass]=alpha*zfs[imass]*betas[imass]/(zas[imass]*betaa);
	
	double *dps=(double*)(ps[imass]),*dr=(double*)r;
	for(int i=0;i<bulk_vol*ndoubles_per_site;i++)
	  {
	    (*dps)=zfs[imass]*(*dr)+alphas[imass]*(*dps);
	    dps++;dr++;
	  }
      }
}

//     -zfs
//     -betas
//     -x
void cgmm2s_update_solutions(basetype **sol,basetype **ps,int *run_flag,double *zfs,double *zas,double *zps,double alpha,double betaa,double betap,double *betas,double *m2,int nmass)
{
  for(int imass=0;imass<nmass;imass++)
    {
      if(run_flag[imass]==1)
	{
	  zfs[imass]=zas[imass]*betap/(betaa*alpha*(1-zas[imass]/zps[imass])+betap*(1-m2[imass]*betaa));
	  betas[imass]=betaa*zfs[imass]/zas[imass];
	  
	  {
	    double *dps=(double*)(ps[imass]),*dsol=(double*)(sol[imass]);
	    for(int i=0;i<bulk_vol*ndoubles_per_site;i++)
	      {
		(*dsol)-=(*dps)*betas[imass];
		dsol++;dps++;
	      }
	  }
	  set_borders_invalid(sol[imass]);
	}
    }
}

void cgmm2s_invert(basetype **sol,quad_su3 **conf,double *m2,int nmass,int niter,double st_res,double st_minres,int st_crit,basetype *source)
{
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
      memset(sol[imass],0,sizeof(basetype)*bulk_vol);
      run_flag[imass]=1;
      set_borders_invalid(sol[imass]);
    }
  
  //     -p=source
  //     -r=source
  //     -calculate Rr=(r,r)
  double rr=cgmm2s_init_vector_and_compute_init_residue(p,r,source);
  
  //writes source norm
  double source_norm=rr;
  if(nissa_inverter_verbosity>=2) master_printf(" Source norm: %lg\n",source_norm);
  if(source_norm==0 || isnan(source_norm)) crash("invalid norm: %lg",source_norm);
  
  //writes initial resiude
  if(nissa_inverter_verbosity>=2)
    {
      master_printf(" cgmms iter 0 rel. residues: ");
      for(int imass=0;imass<nmass;imass++) master_printf("%1.4e  ",1.0);
      master_printf("\n");
    }
  
  //     -betaa=1
  double betaa=1;
  
  //     -ps=source
  //     -zps=zas=1
  //     -alphas=0
  double zps[nmass],zas[nmass],alphas[nmass];
  cgmm2s_init_additional_vectors_and_weights(ps,source,zps,zas,alphas,nmass);
  
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
      if(nrequest!=0) finish_communicating_ev_color_borders(&nrequest,request,p);
      apply_offdiagonal_operator(s,conf,t,p);
      
      //     -pap=(p,s)=(p,Ap)
      double pap=cgmm2s_compute_scalar_product(p,s);
      
      //     calculate betaa=rr/pap=(r,r)/(p,Ap)
      double betap=betaa;
      betaa=-rr/pap;
      
      //     calculate 
      //     -zfs
      //     -betas
      //     -x
      cgmm2s_update_solutions(sol,ps,run_flag,zfs,zas,zps,alpha,betaa,betap,betas,m2,nmass);

      //     calculate
      //     -r'=r+betaa*s=r+beta*Ap
      //     -rfrf=(r',r')
      double rfrf=cgmm2s_compute_new_residue(r,s,betaa);

      //     calculate alpha=rfrf/rr=(r',r')/(r,r)
      alpha=rfrf/rr;
      
      //     calculate p'=r'+alpha*p
      cgmm2s_compute_new_vector(p,r,alpha);
      
      //start the communications of the border
      cgmm2s_start_communicating_borders(&nrequest,request,p);
      
      //     calculate 
      //     -alphas=alpha*zfs*betas/zas*beta
      //     -ps'=r'+alpha*ps
      cgmm2s_scale_additional_vectors(ps,r,alpha,alphas,betaa,betas,zfs,zas,run_flag,nmass);
      
      // shift z and rr
      for(int imass=0;imass<nmass;imass++)
	{
	  zps[imass]=zas[imass];
	  zas[imass]=zfs[imass];
	}
      rr=rfrf;
      iter++;
      
      //     check over residual
      nrun_mass=cgmm2s_check_residue(run_flag,final_res,nrun_mass,rr,zfs,st_crit,st_res,st_minres,iter,nmass,m2,source,conf,s,t,source_norm,sol);
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
	
	w_res=w_res/weight*ndoubles_per_site*bulk_vol*rank_tot;
	max_res*=ndoubles_per_site*bulk_vol*rank_tot;
	
	if(nissa_inverter_verbosity>=1) master_printf(" imass %d, rel residue true=%g approx=%g weighted=%g max=%g\n",imass,res/source_norm,final_res[imass],w_res,max_res);
      }
    }  
  
  if(nissa_inverter_verbosity>=1) master_printf(" Total cgmms iterations: %d\n",iter);
  
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
