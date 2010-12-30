#pragma once

#include "dirac_operator.c"
#include "su3.c"

void inv_Q2_cgmms(spincolor **sol,spincolor *source,spincolor **guess,quad_su3 *conf,double kappac,double *m,int niter,double residue,int nmass)
{
  double zps[nmass],zas[nmass],zfs[nmass],betas[nmass],alphas[nmass];
  double rr,rfrf,pap,alpha;
  double betap,betaa;
  int iter;
  int flag[nmass];

  spincolor *t=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord)); //temporary for internal calculation of DD
  spincolor *s=(spincolor*)malloc(sizeof(spincolor)*(loc_vol));
  spincolor *r=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  spincolor *p=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
  spincolor **ps=(spincolor**)malloc(sizeof(spincolor*)*nmass);
  for(int imass=0;imass<nmass;imass++) ps[imass]=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  
  //sol[*]=0
  for(int imass=0;imass<nmass;imass++)
    {
      memset(sol[imass],0,sizeof(spincolor)*(loc_vol+loc_bord));
      flag[imass]=1;
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
    if(rank_tot>0) MPI_Allreduce(&loc_rr,&rr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    else rr=loc_rr;
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
      apply_Q2(s,p,conf,kappac,0,t);
      
      //     -pap=(p,s)=(p,Ap)
      {
	double loc_pap=0;
	double *dp=(double*)p,*ds=(double*)s;
	for(int i=0;i<loc_vol*3*4*2;i++)
	  {
	    loc_pap+=(*dp)*(*ds);
	    dp++;ds++;
	  }
	if(rank_tot>0) MPI_Allreduce(&loc_pap,&pap,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	else pap=loc_pap;
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
          if(flag[imass]==1)
            {
              zfs[imass]=zas[imass]*zps[imass]*betap/(betaa*alpha*(zps[imass]-zas[imass])+zps[imass]*betap*(1-(m[imass]*m[imass])*betaa));
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
	if(rank_tot>0) MPI_Allreduce(&loc_rfrf,&rfrf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	else rfrf=loc_rfrf;
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

      //     calculate alphas=alpha*zfs*betas/zas*beta
      for(int imass=0;imass<nmass;imass++)
        if(flag[imass]==1)
          {
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
	    
            //     check over residual
            if(rr*zfs[imass]<residue) flag[imass]=0;
            
            //     passes all f into actuals
            zps[imass]=zas[imass];
            zas[imass]=zfs[imass];
          }
      rr=rfrf;
      iter++;

      if(rank==0 && debug) printf("cgmms iter %d residue %g\n",iter,rfrf);
    }
  while(rfrf>residue && iter<niter);

  for(int imass=0;imass<nmass;imass++)  communicate_lx_spincolor_borders(sol[imass]);
  
  for(int imass=0;imass<nmass;imass++) free(ps[imass]);
  free(ps);
  free(s);
  free(p);
  free(r);
  free(t);
}
