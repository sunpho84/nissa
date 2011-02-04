#pragma once

#include "dirac_operator.c"
#include "su3.c"

void inv_Q2_cg_RL(spincolor *sol,spincolor *source,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,int rniter,double residue,int RL)
{
  int riter=0;
  spincolor *s=(spincolor*)malloc(sizeof(spincolor)*(loc_vol));
  spincolor *p=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
  spincolor *r=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  spincolor *t=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord)); //temporary for internal calculation of DD

  if(guess==NULL) memset(sol,0,sizeof(spincolor)*(loc_vol+loc_bord));
  else memcpy(sol,guess,sizeof(spincolor)*(loc_vol+loc_bord));

  //external loop, used if the internal exceed the maximal number of iterations
  double lambda; //(r_(k+1),r_(k+1))
  double source_norm;
  do
    {
      //calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
      double delta;
      {
	apply_Q2_RL(s,sol,conf,kappa,m,t,NULL,NULL,RL);

	double loc_delta=0,loc_source_norm=0;
	double *dsource=(double*)source,*ds=(double*)s,*dp=(double*)p,*dr=(double*)r;
	for(int i=0;i<loc_vol*3*4*2;i++)
	  {
	    double c1=(*dsource)-(*ds);
	    (*dp)=(*dr)=c1;
	    if(riter==0) loc_source_norm+=(*dsource)*(*dsource);
	    loc_delta+=c1*c1;
	    dp++;dr++;ds++;dsource++;
	  }
	MPI_Allreduce(&loc_delta,&delta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	if(riter==0) MPI_Allreduce(&loc_source_norm,&source_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
      

      //main loop
      int iter=0;
      do
	{
	  double omega; //(r_k,r_k)/(p_k*DD*p_k)
	  {
	    double alpha;
	    if(rank_tot>0) communicate_lx_spincolor_borders(p);

	    apply_Q2_RL(s,p,conf,kappa,m,t,NULL,NULL,RL);

	    double loc_alpha=0;
	    complex *cs=(complex*)s,*cp=(complex*)p;
	    for(int i=0;i<loc_vol*3*4;i++)
	      { //real part of the scalar product
		loc_alpha+=((*cs)[0]*(*cp)[0]+(*cs)[1]*(*cp)[1]);

		cs++;cp++;
	      }
	    if(rank_tot>0) MPI_Allreduce(&loc_alpha,&alpha,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	    else alpha=loc_alpha;
	    omega=delta/alpha;
	  }
	  
	  {
	    double loc_lambda=0;
	    double *dsol=(double*)sol,*ds=(double*)s,*dp=(double*)p,*dr=(double*)r;
	    for(int i=0;i<loc_vol*3*4*2;i++)
	      {
		(*dsol)+=omega*(*dp);        //sol_(k+1)=x_k+omega*p_k
		double c1=(*dr)-omega*(*ds); //r_(k+1)=x_k-omega*p_k
		(*dr)=c1;
		loc_lambda+=c1*c1;

		dsol++;ds++;dp++;dr++;
	      }
	    if(rank_tot>0) MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	    else lambda=loc_lambda;
	  }

	  double gammag=lambda/delta; //(r_(k+1),r_(k+1))/(r_k,r_k)
	  delta=lambda;
	  
	  //p_(k+1)=r_(k+1)+gammag*p_k
	  {
	    double *dp=(double*)p,*dr=(double*)r;
	    for(int i=0;i<loc_vol*3*4*2;i++)
	      {
		(*dp)=(*dr)+gammag*(*dp);
		
		dp++;dr++;
	      }
	    }

	  iter++;

	  if(rank==0 && debug && iter%10==0) printf("iter %d residue %g\n",iter,lambda);
	}
      while(lambda>(residue*source_norm) && iter<niter);
      
      //last calculation of residual, in the case iter>niter
      communicate_lx_spincolor_borders(sol);

      apply_Q2_RL(s,sol,conf,kappa,m,t,NULL,NULL,RL);
      {
	double loc_lambda=0;
	double *ds=(double*)s,*dsource=(double*)source;
	for(int i=0;i<loc_vol*3*4*2;i++)
	  {
	    double c1=(*dsource)-(*ds);
	    loc_lambda+=c1*c1;
	    
	    dsource++;ds++;
	  }
	if(rank_tot>0) MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	else lambda=loc_lambda;
      }

      riter++;
    }
  while(lambda>(residue*source_norm) && riter<rniter);
  
  free(s);
  free(p);
  free(r);
  free(t);
}
