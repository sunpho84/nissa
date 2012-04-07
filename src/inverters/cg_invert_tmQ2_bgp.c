#pragma once

#include "base/bgp_instructions.c"

void inv_tmQ2_cg_RL(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,int rniter,double residue,int RL,spincolor *source)
{
  bgp_complex A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32;
  bgp_complex B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32;
  bgp_complex N0,N1,N2,S0,S1,S2;
  bgp_complex R;

  int riter=0;
  spincolor *s=nissa_malloc("s in gc",loc_vol,spincolor);
  spincolor *p=nissa_malloc("p in gc",loc_vol+loc_bord,spincolor);
  spincolor *r=nissa_malloc("r in gc",loc_vol,spincolor);
  spincolor *t=nissa_malloc("t in gc",loc_vol+loc_bord,spincolor); //temporary for internal calculation of DD

  if(guess==NULL) memset(sol,0,sizeof(spincolor)*loc_vol);
  else memcpy(sol,guess,sizeof(spincolor)*loc_vol);
  set_borders_invalid(sol);
  
  //external loop, used if the internal exceed the maximal number of iterations
  double lambda; //(r_(k+1),r_(k+1))
  double source_norm;
  do
    {
      //calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
      double delta;
      {
	apply_tmQ2_RL(s,conf,kappa,m,t,RL,sol);
	
	complex cloc_delta={0,0},cloc_source_norm={0,0};
	
	bgp_color_put_to_zero(N0,N1,N2);
	if(riter==0) bgp_color_put_to_zero(S0,S1,S2);
	nissa_loc_vol_loop(ivol)
	  {
	    bgp_spincolor_load(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,source[ivol]);
	    if(riter==0) bgp_summassign_color_square_spincolor(S0,S1,S2,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);

	    bgp_spincolor_load(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,s[ivol]);
	    bgp_subtassign_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32);
	    bgp_spincolor_save(p[ivol],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	    bgp_spincolor_save(r[ivol],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	    bgp_summassign_color_square_spincolor(N0,N1,N2,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	  }
	set_borders_invalid(p);
	
	bgp_square_norm_color(N0,N1,N2);
	bgp_square_norm_color(S0,S1,S2);
	bgp_complex_save(cloc_delta,N0);
	bgp_complex_save(cloc_source_norm,S0);
		
	MPI_Allreduce(cloc_delta,&delta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(cloc_source_norm,&source_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	if(riter==0) master_printf("Source norm: %lg\n",source_norm);
	master_printf("iter 0 relative residue: 1\n");
      }

      //main loop
      int iter=0;
      do
	{
	  double omega; //(r_k,r_k)/(p_k*DD*p_k)
	  {
	    double alpha;
	    apply_tmQ2_RL(s,conf,kappa,m,t,RL,p);

	    complex cloc_alpha={0,0};

	    bgp_color_put_to_zero(N0,N1,N2);
	    nissa_loc_vol_loop(ivol)
	      { //real part of the scalar product
		bgp_spincolor_load(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,s[ivol]);
		bgp_spincolor_load(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,p[ivol]);
		bgp_summassign_color_realscalarprod_spincolor(N0,N1,N2,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32);
	      }
	    bgp_square_norm_color(N0,N1,N2);
	    bgp_complex_save(cloc_alpha,N0);
	    
	    if(rank_tot>0) MPI_Allreduce(cloc_alpha,&alpha,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	    else alpha=cloc_alpha[0];

	    omega=delta/alpha;
	  }
	  
	  {
	    complex cloc_lambda={0,0};
	    complex comega={omega,0};
	    bgp_complex_load(R,comega);

	    bgp_color_put_to_zero(N0,N1,N2);
	    nissa_loc_vol_loop(ivol)
	      {        //sol_(k+1)=x_k+omega*p_k
		bgp_spincolor_load(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,sol[ivol]);
		bgp_spincolor_load(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,p[ivol]);
		bgp_summassign_spincolor_prod_double(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,omega);
		bgp_spincolor_save(sol[ivol],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);

		//r_(k+1)=x_k-omega*p_k and residue calculation
		bgp_spincolor_load(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,r[ivol]);
		bgp_spincolor_load(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,s[ivol]);
		bgp_subtassign_spincolor_prod_double(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,__creal(R));
		bgp_summassign_color_square_spincolor(N0,N1,N2,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
		bgp_spincolor_save(r[ivol],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	      }
	    bgp_square_norm_color(N0,N1,N2);
	    bgp_complex_save(cloc_lambda,N0);
	    if(rank_tot>0) MPI_Allreduce(cloc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	    else lambda=cloc_lambda[0];
	    set_borders_invalid(sol);
	  }

	  double gammag=lambda/delta; //(r_(k+1),r_(k+1))/(r_k,r_k)
	  delta=lambda;
	  
	  //p_(k+1)=r_(k+1)+gammag*p_k
	  {
	    complex cgammag={gammag,0};
	    bgp_complex_load(R,cgammag);

	    nissa_loc_vol_loop(ivol)
	      {
		bgp_spincolor_load(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,p[ivol]);
		bgp_spincolor_load(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,r[ivol]);
		bgp_summ_spincolor_prod_double(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,__creal(R));
		bgp_spincolor_save(p[ivol],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	      }
	    set_borders_invalid(p);
	    }

	  iter++;

	  if(iter%10==0) master_printf("iter %d relative residue: %g\n",iter,lambda/source_norm);
	}
      while(lambda>(residue*source_norm) && iter<niter);
      
      //last calculation of residual, in the case iter>niter
      apply_tmQ2_RL(s,conf,kappa,m,t,RL,sol);
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
	
	master_printf("\nfinal relative residue (after %d iters): %g where %g was required\n",iter,lambda/source_norm,residue);
      }

      riter++;
    }
  while(lambda>(residue*source_norm) && riter<rniter);
  
  nissa_free(s);
  nissa_free(p);
  nissa_free(r);
  nissa_free(t);
}
