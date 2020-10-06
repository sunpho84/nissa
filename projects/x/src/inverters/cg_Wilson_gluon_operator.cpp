#include <string.h>

#include "../../../../src/nissa.hpp"

#include "../types/types.hpp"
#include "../operators/Wilson_gluon_Klein_Gordon_operator.hpp"

void inv_Wilson_gluon_Klein_Gordon_operator(spin1field *sol,spin1field *guess,gluon_info gl,int niter,int rniter,double residue,spin1field *source)
{
  int riter=0;
  spin1field *s=nissa_malloc("s",loc_vol,spin1field);
  spin1field *p=nissa_malloc("p",loc_vol+bord_vol,spin1field);
  spin1field *r=nissa_malloc("r",loc_vol,spin1field);

  if(guess==NULL) memset(sol,0,sizeof(spin1field)*(loc_vol+bord_vol));
  else memcpy(sol,guess,sizeof(spin1field)*(loc_vol+bord_vol));
  set_borders_invalid(sol);  
  
  //external loop, used if the internal exceed the maximal number of iterations
  double lambda; //(r_(k+1),r_(k+1))
  double source_norm;
  do
    {
      //calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
      double delta;
      {
	apply_Wilson_gluon_x_Klein_Gordon_operator(s,sol,gl);
	
	double loc_delta=0,loc_source_norm=0;
	double *dsource=(double*)source,*ds=(double*)s,*dp=(double*)p,*dr=(double*)r;
	for(int i=0;i<loc_vol*4*2;i++)
	  {
	    double c1=(*dsource)-(*ds);
	    (*dp)=(*dr)=c1;
	    loc_source_norm+=(*dsource)*(*dsource);
	    loc_delta+=c1*c1;
	    dp++;dr++;ds++;dsource++;
	  }
	MPI_Allreduce(&loc_delta,&delta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&loc_source_norm,&source_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	if(riter==0) verbosity_lv2_master_printf("Source norm: %lg\n",source_norm);
	verbosity_lv2_master_printf("iter 0 relative residue: %lg\n",delta/source_norm);
	set_borders_invalid(p);
      }

      //main loop
      int iter=0;
      do
	{
	  double omega; //(r_k,r_k)/(p_k*DD*p_k)
	  {
	    double alpha;

	    apply_Wilson_gluon_x_Klein_Gordon_operator(s,p,gl);
	    
	    double loc_alpha=0;
	    complex *cs=(complex*)s,*cp=(complex*)p;
	    for(int i=0;i<loc_vol*4;i++)
	      { //real part of the scalar product
		loc_alpha+=((*cs)[0]*(*cp)[0]+(*cs)[1]*(*cp)[1]);

		cs++;cp++;
	      }
	    if(nranks>0) MPI_Allreduce(&loc_alpha,&alpha,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	    else alpha=loc_alpha;
	    omega=delta/alpha;
	  }
	  
	  {
	    double loc_lambda=0;
	    double *dsol=(double*)sol,*ds=(double*)s,*dp=(double*)p,*dr=(double*)r;
	    for(int i=0;i<loc_vol*4*2;i++)
	      {
		(*dsol)+=omega*(*dp);        //sol_(k+1)=x_k+omega*p_k
		double c1=(*dr)-omega*(*ds); //r_(k+1)=x_k-omega*p_k
		(*dr)=c1;
		loc_lambda+=c1*c1;

		dsol++;ds++;dp++;dr++;
	      }
	    if(nranks>0) MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	    else lambda=loc_lambda;
	    set_borders_invalid(sol);
	  }

	  double gammag=lambda/delta; //(r_(k+1),r_(k+1))/(r_k,r_k)
	  delta=lambda;
	  
	  //p_(k+1)=r_(k+1)+gammag*p_k
	  {
	    double *dp=(double*)p,*dr=(double*)r;
	    for(int i=0;i<loc_vol*4*2;i++)
	      {
		(*dp)=(*dr)+gammag*(*dp);
		
		dp++;dr++;
	      }
	    set_borders_invalid(p);
	    }

	  iter++;

	  if(iter%10==0) verbosity_lv2_master_printf("iter %d relative residue: %lg\n",iter,lambda/source_norm);
	}
      while(lambda>(residue*source_norm) && iter<niter);
      
      //last calculation of residual, in the case iter>niter
      
      apply_Wilson_gluon_x_Klein_Gordon_operator(s,sol,gl);
      {
	double loc_lambda=0;
	double *ds=(double*)s,*dsource=(double*)source;
	for(int i=0;i<loc_vol*4*2;i++)
	  {
	    double c1=(*dsource)-(*ds);
	    loc_lambda+=c1*c1;
	    
	    dsource++;ds++;
	  }
	if(nranks>0) MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	else lambda=loc_lambda;
	
	verbosity_lv1_master_printf("\nfinal relative residue (after %d iters): %lg where %lg was required\n",iter,lambda/source_norm,residue);
      }

      riter++;
    }
  while(lambda>(residue*source_norm) && riter<rniter);
  
  nissa_free(s);
  nissa_free(p);
  nissa_free(r);
}
