#include <string.h>

#include "../../new_types/new_types_definitions.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../base/routines.h"
#include "../../geometry/geometry_mix.h"
#include "../../dirac_operators/dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr.h"
#include "cg_128_invert_tmDeoimpr.h"

//Refers to the doc: "doc/eo_inverter.lyx" for explenations

//invert Koo defined in equation (7)
void inv_tmDkern_eoprec_square_eos(spincolor *sol,spincolor *guess,quad_su3 **conf,double kappa,double mu,int nitermax,double residue,spincolor *source)
{
  int niter=nitermax;
  int riter=0;
  int rniter=5;
  spincolor *p=nissa_malloc("p",loc_volh+bord_volh,spincolor);
  spincolor *r=nissa_malloc("r",loc_volh,spincolor);
  spincolor *s=nissa_malloc("s",loc_volh,spincolor);
  spincolor *temp1=nissa_malloc("temp1",loc_volh+bord_volh,spincolor);
  spincolor *temp2=nissa_malloc("temp2",loc_volh+bord_volh,spincolor);
  
  ///////////////// prepare the internal source /////////////////
  
  if(guess==NULL) memset(sol,0,sizeof(spincolor)*loc_volh);
  else memcpy(sol,guess,sizeof(spincolor)*loc_volh);
  set_borders_invalid(sol);
  
  //external loop, used if the internal exceed the maximal number of iterations
  double lambda; //(r_(k+1),r_(k+1))
  double source_norm;
  do
    {
      //calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
      double delta;
      {
	tmDkern_eoprec_square_eos(s,temp1,temp2,conf,kappa,mu,sol);
	
        double loc_delta=0,loc_source_norm=0;
	nissa_loc_volh_loop(ivol)
	  for(int id=0;id<4;id++)
	    for(int ic=0;ic<3;ic++)
	      for(int ri=0;ri<2;ri++)
		{
		  double c1=source[ivol][id][ic][ri]-s[ivol][id][ic][ri];
		  p[ivol][id][ic][ri]=r[ivol][id][ic][ri]=c1;
		  if(riter==0) loc_source_norm+=source[ivol][id][ic][ri]*source[ivol][id][ic][ri];
		  loc_delta+=c1*c1;
		}	
	set_borders_invalid(p);
	
        MPI_Allreduce(&loc_delta,&delta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	if(riter==0)
	  {
	    MPI_Allreduce(&loc_source_norm,&source_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	    master_printf("\nSource norm: %lg\n",source_norm);
	    master_printf("iter 0 relative residue: %lg\n",delta/source_norm);
	  }
      }
      
      //main loop
      int iter=0;
      do
	{
	  double omega; //(r_k,r_k)/(p_k*DD*p_k)
	  double alpha;
	  
	  tmDkern_eoprec_square_eos(s,temp1,temp2,conf,kappa,mu,p);
	  
	  double loc_alpha=0; //real part of the scalar product
	  nissa_loc_volh_loop(ivol)
	    for(int id=0;id<4;id++)
	      for(int ic=0;ic<3;ic++)
		for(int ri=0;ri<2;ri++)
		  loc_alpha+=s[ivol][id][ic][ri]*p[ivol][id][ic][ri];
	  MPI_Allreduce(&loc_alpha,&alpha,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  omega=delta/alpha;
	  
	  double loc_lambda=0;
	  nissa_loc_volh_loop(ivol)
	    for(int id=0;id<4;id++)
	      for(int ic=0;ic<3;ic++)
		for(int ri=0;ri<2;ri++)
		  {
		    sol[ivol][id][ic][ri]+=omega*p[ivol][id][ic][ri];    //sol_(k+1)=x_k+omega*p_k
		    double c1=r[ivol][id][ic][ri]-omega*s[ivol][id][ic][ri];//r_(k+1)=x_k-omega*pk
		    r[ivol][id][ic][ri]=c1;
		    loc_lambda+=c1*c1;
		  }
	  MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  set_borders_invalid(sol);
	  
	  double gammag=lambda/delta; //(r_(k+1),r_(k+1))/(r_k,r_k)
	  delta=lambda;
	  
	  //p_(k+1)=r_(k+1)+gammag*p_k
	  nissa_loc_volh_loop(ivol)
	    for(int id=0;id<4;id++)
	      for(int ic=0;ic<3;ic++)
		for(int ri=0;ri<2;ri++)
		  p[ivol][id][ic][ri]=r[ivol][id][ic][ri]+gammag*p[ivol][id][ic][ri];
	  set_borders_invalid(p);
	  
	  iter++;

          if(iter%10==0) master_printf("iter %d relative residue: %lg\n",iter,lambda/source_norm);
	}
      while(lambda>(residue*source_norm) && iter<niter);
      
      //last calculation of residual, in the case iter>niter
      tmDkern_eoprec_square_eos(s,temp1,temp2,conf,kappa,mu,sol);
      {
        double loc_lambda=0;
	nissa_loc_volh_loop(ivol)
	  for(int id=0;id<4;id++)
	    for(int ic=0;ic<3;ic++)
	      for(int ri=0;ri<2;ri++)
		{
		  double c1=source[ivol][id][ic][ri]-s[ivol][id][ic][ri];
		  loc_lambda+=c1*c1;
		}
	if(nissa_nranks>0) MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	else lambda=loc_lambda;
      }
      master_printf("\nfinal relative residue (after %d iters): %lg where %lg was required\n",iter,lambda/source_norm,residue);

      riter++;
    }
  while(lambda>(residue*source_norm) && riter<rniter);

  nissa_free(s);
  nissa_free(p);
  nissa_free(r);
  nissa_free(temp1);
  nissa_free(temp2);
}

//hack for template, which require rniter
void inv_tmDkern_eoprec_square_eos(spincolor *sol,spincolor *guess,quad_su3 **conf,double kappa,double mu,int nitermax,int rniter,double residue,spincolor *source)
{inv_tmDkern_eoprec_square_eos(sol,guess,conf,kappa,mu,nitermax,residue,source);}

//Invert twisted mass operator using e/o preconditioning.
void inv_tmD_cg_eoprec_eos(spincolor *solution_lx,spincolor *guess_Koo,quad_su3 *conf_lx,double kappa,double mu,int nitermax,double residue,spincolor *source_lx)
{
  //prepare the e/o split version of the source
  spincolor *source_eos[2];
  source_eos[0]=nissa_malloc("source_eos0",loc_volh+bord_volh,spincolor);
  source_eos[1]=nissa_malloc("source_eos1",loc_volh+bord_volh,spincolor);
  split_lx_spincolor_into_eo_parts(source_eos,source_lx);
  
  //prepare the e/o split version of the solution
  spincolor *solution_eos[2];
  solution_eos[0]=nissa_malloc("solution_eos_0",loc_volh+bord_volh,spincolor);
  solution_eos[1]=nissa_malloc("solution_eos_0",loc_volh+bord_volh,spincolor);
  
  //prepare the e/o split version of the conf
  quad_su3 *conf_eos[2];
  conf_eos[0]=nissa_malloc("conf_eos_0",loc_volh+bord_volh,quad_su3);
  conf_eos[1]=nissa_malloc("conf_eos_1",loc_volh+bord_volh,quad_su3);
  split_lx_conf_into_eo_parts(conf_eos,conf_lx);
  
  ///////////////////////////////////// invert with e/o improvement ///////////////////////////////////
  
  spincolor *varphi=nissa_malloc("varphi",loc_volh+bord_volh,spincolor);
  
  //Equation (8.a)
  spincolor *temp=nissa_malloc("temp",loc_volh+bord_volh,spincolor);
  inv_tmDee_or_oo_eos(temp,kappa,mu,source_eos[EVN]);
  
  //Equation (8.b)
  tmn2Doe_eos(varphi,conf_eos,temp);
  nissa_loc_volh_loop(ivol)
    for(int id=0;id<2;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  { //gamma5 is explicitely wrote
	    varphi[ivol][id  ][ic][ri]=+source_eos[ODD][ivol][id  ][ic][ri]+varphi[ivol][id  ][ic][ri]*0.5;
	    varphi[ivol][id+2][ic][ri]=-source_eos[ODD][ivol][id+2][ic][ri]-varphi[ivol][id+2][ic][ri]*0.5;
	  }
  set_borders_invalid(varphi);
  
  //Equation (9) using solution_eos[EVN] as temporary vector
  if(nissa_use_128_bit_precision) inv_tmDkern_eoprec_square_eos_128(temp,guess_Koo,conf_eos,kappa,mu,nitermax,residue,varphi);
  else inv_tmDkern_eoprec_square_eos(temp,guess_Koo,conf_eos,kappa,mu,nitermax,residue,varphi);
  tmDkern_eoprec_eos(solution_eos[ODD],solution_eos[EVN],conf_eos,kappa,-mu,temp);
  if(guess_Koo!=NULL) memcpy(guess_Koo,temp,sizeof(spincolor)*loc_volh); //if a guess was passed, return new one
  nissa_free(temp);

  //Equation (10)
  tmn2Deo_eos(varphi,conf_eos,solution_eos[ODD]);
  nissa_loc_volh_loop(ivol)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  varphi[ivol][id][ic][ri]=source_eos[EVN][ivol][id][ic][ri]+varphi[ivol][id][ic][ri]*0.5;
  set_borders_invalid(varphi);
  inv_tmDee_or_oo_eos(solution_eos[EVN],kappa,mu,varphi);

  nissa_free(varphi);

  /////////////////////////// paste the e/o parts of the solution together and free ///////////////////
  
  paste_eo_parts_into_lx_spincolor(solution_lx,solution_eos);
  
  for(int par=0;par<2;par++)
    {
      nissa_free(conf_eos[par]);
      nissa_free(source_eos[par]);
      nissa_free(solution_eos[par]);
    }
}
