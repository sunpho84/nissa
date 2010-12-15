#pragma once

void inv_DDdag_cg(spincolor *sol,spincolor *source,quad_su3 *conf,double kappac,double m,int niter,int rniter)
{
  int riter=0;
  spincolor *s=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  spincolor *p=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  spincolor *r=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  spincolor *t=(spincolor*)malloc(sizeof(spincolor)*loc_vol); //temporary for internal calculation of DD

  put_to_zero_locvol_spincolor(sol);
  
  //external loop, used if the internal exceed the maximal number of iterations
  do
    {
      //calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
      double delta;
      {
	apply_DDdag(s,sol,kappac,m,t);
	double loc_delta=0;
	double *dsource=(double*)source,*ds=(double*)s,*dp=(double*)p,*dr=(double*)r;
	for(int i=0;i<loc_vol*3*2;i++)
	  {
	    double c1=(*dsource)-(*ds);
	    (*dp)=(*dr)=c1;
	    loc_delta+=c1*c1;
	    
	    dp++;dr++;ds++;dsource++;
	  }
	MPI_Allreduce(&loc_delta,&glb_delta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }

      //main loop
      int iter=0;
      do
	{

	  double omega; //(r_k,r_k)/(p_k*DD*p_k)
	  {
	    double alpha;
	    apply_DDdag(s,p,kappac,m,t);
	    double loc_alpha=0;
	    complex *cs=(complex*)s,*cp=(complex*)p;
	    for(int i=0;i<loc_vol*3;i++)
	      { //real part of the scalar product
		loc_alpha+=((*cs)[0]*(*cp)[0]+(*cs)[1]*(*cp)[1]);

		cs++;cp++;
	      }
	    MPI_Allreduce(&loc_alpha,&alpha,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	    omega=delta/alpha;
	  }
	  
	  double lambda; //(r_(k+1),r_(k+1))
	  {
	    double loc_lambda=0;
	    double *dsol=(double*)sol,*ds=(double*)s,*dp=(double*)p,*dr=(double*)r;
	    for(int i=0;i<loc_vol*3*2;i++)
	      {
		(*dsol)+=omega*(*dp);        //sol_(k+1)=x_k+omega*p_k
		double c1=(*dr)-omega*(*ds); //r_(k+1)=x_k-omega*p_k
		(*dr)=c1;
		loc_lambda+=c1*c1;

		dsol++;ds++;dp++;dr++;
	      }
	    MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  }

	  double gammag=lambda/delta; //(r_(k+1),r_(k+1))/(r_k,r_k)
	  delta=lambda;
	  
	  //p_(k+1)=r_(k+1)+gammag*p_k
	  {
	    double *dp=(double*)p,*dr=(double*)r;
	    for(int i=0;i<loc_vol*3*2;i++)
	      {
		(*dp)=(*r)+gammag*(*p);
		
		dp++;dr++;
	      }
	    }

	  iter++;

	}
      while(lambda>residue && iter<niter);
      
      //last calculation of residual, in the case iter>niter
      double lambda;
      apply_DDdag(s,sol,kappac,m,t);
      {
	double loc_lambda=0;
	double *ds=(double*)s,*dsource=(double*)source;
	for(int i=0;i<loc_vol*3*2;i++)
	  {
	    double c1=(*dsource)-(*ds);
	    loc_lambda+=c1*c1;
	    
	    dsource++;ds++;
	  }
	MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }

      riter++;
    }
  while(lambda>residue && riter<rniter);
  
  free(s);
  free(p);
  free(r);
  free(t);
}
