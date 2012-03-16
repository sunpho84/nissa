#pragma once

void inv_stD2ee_cgmm2s(color **sol,color *source,quad_su3 **conf,double *m2,int nmass,int niter,double st_res,double st_minres,int st_crit)
{
  double zps[nmass],zas[nmass],zfs[nmass],betas[nmass],alphas[nmass];
  double rr,rfrf,pap,alpha;
  double betap,betaa;
  int iter;
  int run_flag[nmass],nrun_mass=nmass;
  double final_res[nmass];
  double source_norm;
  
  color *t=nissa_malloc("DD_temp",loc_volh+loc_bordh,color);
  color *s=nissa_malloc("s",loc_volh,color);
  color *r=nissa_malloc("r",loc_volh,color);
  color *p=nissa_malloc("p",loc_volh+loc_bordh,color);
  color *ps[nmass];
  for(int imass=0;imass<nmass;imass++) ps[imass]=nissa_malloc("ps",loc_volh,color);
  
  communicate_eo_quad_su3_borders(conf[0],conf[1]);
  
  //sol[*]=0
  for(int imass=0;imass<nmass;imass++)
    {
      memset(sol[imass],0,sizeof(color)*(loc_volh+loc_bordh));
      run_flag[imass]=1;
    }
  
  //     -p=source
  //     -r=source
  //     -calculate Rr=(r,r)
  {
    double loc_rr=0;
    double *dsource=(double*)source,*dp=(double*)p,*dr=(double*)r;
    for(int i=0;i<loc_volh*3*2;i++)
      {
	(*dp)=(*dsource);
	(*dr)=(*dsource);
	loc_rr+=(*dr)*(*dr);
	
	dsource++;dp++;dr++;
      }
    MPI_Allreduce(&loc_rr,&rr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    source_norm=rr;
    master_printf(" Source norm: %lg\n",source_norm);
    
    master_printf(" cgmms iter 0 rel. residues: ");
    for(int imass=0;imass<nmass;imass++) master_printf("%1.4e  ",1.0);
    master_printf("\n");
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
	for(int i=0;i<loc_volh*3*2;i++)
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
      communicate_ev_color_borders(p);
      apply_st2Doe(t,conf,p);
      communicate_od_color_borders(t);
      apply_st2Deo(s,conf,t);
      
      //     -pap=(p,s)=(p,Ap)
      {
	double loc_pap=0;
	double *dp=(double*)p,*ds=(double*)s;
	int n=0;
	for(int i=0;i<loc_volh*3*2;i++)
	  {
	    (*ds)*=0.25;
	    loc_pap+=(*dp)*(*ds);
	    dp++;ds++;
	    n++;
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
              zfs[imass]=zas[imass]*betap/(betaa*alpha*(1-zas[imass]/zps[imass])+betap*(1-m2[imass]*betaa));
              betas[imass]=betaa*zfs[imass]/zas[imass];
	      
	      {
		double *dps=(double*)(ps[imass]),*dsol=(double*)(sol[imass]);
		for(int i=0;i<loc_volh*3*2;i++)
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
	for(int i=0;i<loc_volh*3*2;i++)
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
	for(int i=0;i<loc_volh*3*2;i++)
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
	      for(int i=0;i<loc_volh*3*2;i++)
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
      nrun_mass=check_cgmm2s_residue_stD2ee(run_flag,final_res,nrun_mass,rr,zfs,st_crit,st_res,st_minres,iter,sol,nmass,m2,source,conf,s,t,source_norm);
    }
  while(nrun_mass>0 && iter<niter);
  
  for(int imass=0;imass<nmass;imass++) communicate_ev_color_borders(sol[imass]);
  
  //print the final true residue
  
  for(int imass=0;imass<nmass;imass++)
    {
      double res,w_res,weight,max_res;
      apply_stD2ee(s,conf,t,sqrt(m2[imass]),sol[imass]);
      {
	double loc_res=0;
	double locw_res=0;
	double locmax_res=0;
	double loc_weight=0;
	
	complex *ds=(complex*)s,*dsour=(complex*)source,*dsol=(complex*)(sol[imass]);
	for(int i=0;i<loc_volh*3;i++)
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
	
	master_printf(" imass %d, rel residue true=%g approx=%g weighted=%g max=%g\n",imass,res/source_norm,final_res[imass],w_res,max_res);
      }
    }  
  
  master_printf(" Total cgmms iterations: %d\n",iter);
  
  for(int imass=0;imass<nmass;imass++) nissa_free(ps[imass]);
  nissa_free(s);
  nissa_free(p);
  nissa_free(r);
  nissa_free(t);
}
