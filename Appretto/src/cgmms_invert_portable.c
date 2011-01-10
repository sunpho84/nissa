#pragma once

#include "dirac_operator.c"
#include "su3.c"

double calculate_weighted_residue(spincolor *source,spincolor *sol,quad_su3 *conf,double kappac,double m,spincolor *s,spincolor *t,int dinf,redspincolor *tin,redspincolor *tout)
{
  apply_Q2(s,sol,conf,kappac,m,t,tout,tin);
  
  double loc_weighted_residue=0,tot_weighted_residue;
  double loc_weight=0,tot_weight;
  double *ds=(double*)s,*dsource=(double*)source,*dsol=(double*)sol;
  
  for(int i=0;i<loc_vol*3*4;i++)
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

  if(rank_tot>0)
    if(dinf==2)
      {
	MPI_Allreduce(&loc_weighted_residue,&tot_weighted_residue,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&loc_weight,&tot_weight,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
    else MPI_Allreduce(&loc_weighted_residue,&tot_weighted_residue,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  else tot_weighted_residue=loc_weighted_residue;
  
  if(dinf==2) tot_weighted_residue/=loc_vol*tot_weight;
  
  return tot_weighted_residue;
}

void inv_Q2_cgmms(spincolor **sol,spincolor *source,spincolor **guess,quad_su3 *conf,double kappac,double *m,int nmass,int niter,double stopping_residue,int stopping_criterion)
{
  double zps[nmass],zas[nmass],zfs[nmass],betas[nmass],alphas[nmass];
  double rr,rfrf,pap,alpha;
  double betap,betaa;
  int iter;
  int flag[nmass],stop=0,running_mass=nmass;

  redspincolor *tout=(redspincolor*)malloc(sizeof(redspincolor)*loc_bord);
  redspincolor *tin=(redspincolor*)malloc(sizeof(redspincolor)*loc_bord);

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

    if(stopping_criterion==sc_standard||stopping_criterion==sc_differentiate) stopping_residue*=rr;
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
      if(rank_tot>0) communicate_lx_redspincolor_borders(p,tout,tin);
      apply_Q2(s,p,conf,kappac,m[0],t,tout,tin);
      
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

      for(int imass=0;imass<nmass;imass++)
        if(flag[imass]==1)
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
      double residue;
      if(stopping_criterion==sc_standard)
	{
	  residue=rfrf;
	  if(residue<stopping_residue) stop=1;

	  if(rank==0 && debug) printf("cgmms iter %d residue %g\n",iter,residue);
     	}
      else	      //different stopping criterion for each mass
	{
	  if((stopping_criterion==sc_differentiate||iter%10==0) && rank==0 && debug)
	    printf("cgmms iter %d residue(s):\t",iter);
	for(int imass=nmass-1;imass>=0;imass--)
	  if(flag[imass])
	    {
	      if(stopping_criterion==sc_differentiate)
		{
		  residue=rr*zfs[imass];
		  if(rank==0 && debug) printf("%g\t",residue);
		}
	      else
		if(iter%10==0)
		  {  //locally weighted norm
		    communicate_lx_redspincolor_borders(sol[imass],tout,tin);
		    if(stopping_criterion==sc_weighted_norm2)
		      residue=calculate_weighted_residue(source,sol[imass],conf,kappac,m[imass],s,t,2,tout,tin);
		    else
		      residue=calculate_weighted_residue(source,sol[imass],conf,kappac,m[imass],s,t,-1,tout,tin);
		  	      
		    if(residue<stopping_residue)
		      {
			running_mass--;
			flag[imass]=0;
		      }
		    if(rank==0 && debug) printf("%g\t",residue);
		  }
	    }
	if((stopping_criterion==sc_differentiate||iter%10==0) && rank==0 && debug)
	  printf("\n");
	}
      if(running_mass==0) stop=1;
    }
  while(stop==0 && iter<niter);

  for(int imass=0;imass<nmass;imass++)  communicate_lx_spincolor_borders(sol[imass]);
  
  for(int imass=0;imass<nmass;imass++) free(ps[imass]);
  free(ps);
  free(s);
  free(p);
  free(r);
  free(t);

  free(tout);
  free(tin);
}
