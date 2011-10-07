#pragma once

void inv_Q2_cgmms_RL(spincolor **sol,spincolor *source,spincolor **guess,quad_su3 *conf,double kappa,double *m,int nmass,int niter,double st_res,double st_minres,int st_crit,int RL)
{
  double zps[nmass],zas[nmass],zfs[nmass],betas[nmass],alphas[nmass];
  double rr,rfrf,pap,alpha;
  double betap,betaa;
  int iter;
  int run_flag[nmass],nrun_mass=nmass;
  double final_res[nmass];

  spincolor *t=appretto_malloc("DD_temp",loc_vol+loc_bord,spincolor);
  spincolor *s=appretto_malloc("s",loc_vol,spincolor);
  spincolor *r=appretto_malloc("r",loc_vol,spincolor);
  spincolor *p=appretto_malloc("p",loc_vol+loc_bord,spincolor);
  spincolor *ps[nmass];
  for(int imass=0;imass<nmass;imass++) ps[imass]=appretto_malloc("ps",loc_vol,spincolor);

  //sol[*]=0
  for(int imass=0;imass<nmass;imass++)
    {
      memset(sol[imass],0,sizeof(spincolor)*(loc_vol+loc_bord));
      run_flag[imass]=1;
    }
  
  //     -p=source
  //     -r=source
  //     -calculate Rr=(r,r)
  {
    double loc_rr=0;
    double *dsource=(double*)source,*dp=(double*)p,*dr=(double*)r;
    for(int i=0;i<loc_vol*3*4*2;i++)
      {
	(*dp)=(*dsource);
	(*dr)=(*dsource);
	loc_rr+=(*dr)*(*dr);
	
	dsource++;dp++;dr++;
      }
    MPI_Allreduce(&loc_rr,&rr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    if(st_crit==sc_standard||st_crit==sc_unilevel) st_res*=rr;
    
    printf("cgmms iter 0 residues: ");
    for(int imass=0;imass<nmass;imass++) printf("%1.4e  ",rr);
    printf("\n");
  }

  //     -betaa=1
  betaa=1;
  
  //     -ps=source
  //     -zps=zas=1
  //     -alphas=0
  {
    for(int imass=0;imass<nmass;imass++)
      {
	double *dps=(double*)(ps[imass]),*dsource=(double*)source;
	for(int i=0;i<loc_vol*3*4*2;i++)
	  {
	    (*dps)=(*dsource);
	    dps++;dsource++;
	  }

	zps[imass]=zas[imass]=1;
	alphas[imass]=0;
      }
  }

  //     -alpha=0
  alpha=0;
      
  iter=0;
  do
    {
      //     -s=Ap
      if(rank_tot>0) communicate_lx_spincolor_borders(p);
      apply_Q2_RL(s,p,conf,kappa,m[0],t,NULL,NULL,RL);
      
      //     -pap=(p,s)=(p,Ap)
      {
	double loc_pap=0;
	double *dp=(double*)p,*ds=(double*)s;
	for(int i=0;i<loc_vol*3*4*2;i++)
	  {
	    loc_pap+=(*dp)*(*ds);
	    dp++;ds++;
	  }
	MPI_Allreduce(&loc_pap,&pap,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }

      //     calculate betaa=rr/pap=(r,r)/(p,Ap)
      betap=betaa;
      betaa=-rr/pap;
      
      //     calculate 
      //     -zfs
      //     -betas
      //     -x
      for(int imass=0;imass<nmass;imass++)
        {
          if(run_flag[imass]==1)
            {
              zfs[imass]=zas[imass]*zps[imass]*betap/(betaa*alpha*(zps[imass]-zas[imass])+zps[imass]*betap*(1-(m[imass]-m[0])*(m[imass]+m[0])*betaa));
              betas[imass]=betaa*zfs[imass]/zas[imass];
	      
	      {
		double *dps=(double*)(ps[imass]),*dsol=(double*)(sol[imass]);
		for(int i=0;i<loc_vol*3*4*2;i++)
		  {
		    (*dsol)-=(*dps)*betas[imass];
		    dsol++;dps++;
		  }
	      }
            }
        }
      
      //     calculate
      //     -r'=r+betaa*s=r+beta*Ap
      //     -rfrf=(r',r')
      {
	double loc_rfrf=0;

	double *dr=(double*)r,*ds=(double*)s;	
	for(int i=0;i<loc_vol*3*4*2;i++)
	  {
            (*dr)+=(*ds)*betaa;
            loc_rfrf+=(*dr)*(*dr);
	    dr++;ds++;
          }
	MPI_Allreduce(&loc_rfrf,&rfrf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }

      //     calculate alpha=rfrf/rr=(r',r')/(r,r)
      alpha=rfrf/rr;
      
      //     calculate p'=r'+alpha*p
      {
	double *dp=(double*)p,*dr=(double*)r;		
	for(int i=0;i<loc_vol*3*4*2;i++)
          {
            (*dp)=(*dr)+alpha*(*dp);
	    dp++;dr++;
          }
      }

      for(int imass=0;imass<nmass;imass++)
        if(run_flag[imass]==1)
          { //     calculate alphas=alpha*zfs*betas/zas*beta
            alphas[imass]=alpha*zfs[imass]*betas[imass]/(zas[imass]*betaa);
            
            //     calculate ps'=r'+alpha*p
	    {
	      double *dps=(double*)(ps[imass]),*dr=(double*)r,*ds=(double*)s;
	      for(int i=0;i<loc_vol*3*4*2;i++)
                {
                  (*dps)=zfs[imass]*(*dr)+alphas[imass]*(*dps);
		  dps++;ds++;dr++;
                }
	    }
	    
            //     passes all f into actuals
            zps[imass]=zas[imass];
            zas[imass]=zfs[imass];
          }
      rr=rfrf;
      
      iter++;
      
      //     check over residual
      nrun_mass=check_cgmms_residue_RL(run_flag,final_res,nrun_mass,rr,zfs,st_crit,st_res,st_minres,iter,sol,nmass,m,source,conf,kappa,s,t,RL);
      
    }
  while(nrun_mass>0 && iter<niter);
  
  for(int imass=0;imass<nmass;imass++) communicate_lx_spincolor_borders(sol[imass]);
  
  //print the final true residue
  
  for(int imass=0;imass<nmass;imass++)
    {
      double res,w_res,weight,max_res;
      apply_Q2_RL(s,sol[imass],conf,kappa,m[imass],t,NULL,NULL,RL);
      {
	double loc_res=0;
	double locw_res=0;
	double locmax_res=0;
	double loc_weight=0;

	complex *ds=(complex*)s,*dsour=(complex*)source,*dsol=(complex*)(sol[imass]);
	for(int i=0;i<loc_vol*3*4;i++)
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

	MPI_Reduce(&loc_res,&res,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&locw_res,&w_res,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&loc_weight,&weight,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&locmax_res,&max_res,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

	w_res=w_res/weight*12*glb_vol;
	max_res*=12*glb_vol;
	
	master_printf("imass %d, residue true=%g approx=%g weighted=%g max=%g\n",imass,res,final_res[imass],w_res,max_res);
      }
    }  

  for(int imass=0;imass<nmass;imass++) appretto_free(ps[imass]);
  appretto_free(s);
  appretto_free(p);
  appretto_free(r);
  appretto_free(t);
}
