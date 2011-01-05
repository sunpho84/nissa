#pragma once

#include "dirac_operator.c"
#include "su3.c"

void inv_Q2_cgmms(spincolor **sol,spincolor *source,spincolor **guess,quad_su3 *conf,double kappac,double *m,int nmass,int niter,double stopping_residue,int stopping_criterion)
{
  static double _Complex A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32;
  static double _Complex B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32;
  static double _Complex N0,N1,N2;
  static double _Complex R;

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
    complex cloc_rr={0,0};
    bgp_color_put_to_zero(N0,N1,N2);

    for(int i=0;i<loc_vol;i++)
      {
	bgp_load_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,source[i]);
	bgp_save_spincolor(p[i],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	bgp_save_spincolor(r[i],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	bgp_summassign_color_square_spincolor(N0,N1,N2,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
      }
    bgp_square_norm_color(N0,N1,N2);
    bgp_save_complex(cloc_rr,N0);
        
    if(rank_tot>0) MPI_Allreduce(cloc_rr,&rr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    else rr=cloc_rr[0];
  }

  //     -betaa=1
  betaa=1;
  
  //     -ps=source
  {
    for(int i=0;i<loc_vol;i++)
      {
	bgp_load_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,source[i]);
	for(int imass=0;imass<nmass;imass++)
	  bgp_save_spincolor(ps[imass][i],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
      }
  }

  //     -zps=zas=1
  //     -alphas=0
  {
    for(int imass=0;imass<nmass;imass++)
      {
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
      apply_Q2(s,p,conf,kappac,m[0],t);
      
      //     -pap=(p,s)=(p,Ap)
      {
	complex cloc_pap={0,0};
        bgp_color_put_to_zero(N0,N1,N2);

	for(int i=0;i<loc_vol;i++)
	  {
	    bgp_load_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,s[i]);
            bgp_load_spincolor(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,p[i]);
	    bgp_summassign_color_realscalarprod_spincolor(N0,N1,N2,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32);
	  }
	bgp_square_norm_color(N0,N1,N2);
	bgp_save_complex(cloc_pap,N0);

	if(rank_tot>0) MPI_Allreduce(cloc_pap,&pap,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	else pap=cloc_pap[0];
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
		for(int i=0;i<loc_vol;i++)
		  {
		    bgp_load_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,sol[imass][i]);
		    bgp_load_spincolor(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,ps[imass][i]);
		    bgp_subtassign_spincolor_prod_real(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,betas[imass]);
		    bgp_save_spincolor(sol[imass][i],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
		  }
	      }
            }
        }
      
      //     calculate
      //     -r'=r+betaa*s=r+beta*Ap
      //     -rfrf=(r',r')
      {
	complex cloc_rfrf={0,0};
        bgp_color_put_to_zero(N0,N1,N2);

	for(int i=0;i<loc_vol;i++)
	  {
	    bgp_load_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,r[i]);
	    bgp_load_spincolor(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,s[i]);
	    bgp_summassign_spincolor_prod_real(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,betaa);
	    bgp_summassign_color_square_spincolor(N0,N1,N2,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	    bgp_save_spincolor(r[i],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
          }
	bgp_square_norm_color(N0,N1,N2);
	bgp_save_complex(cloc_rfrf,N0);

	if(rank_tot>0) MPI_Allreduce(cloc_rfrf,&rfrf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	else rfrf=cloc_rfrf[0];
      }

      //     calculate alpha=rfrf/rr=(r',r')/(r,r)
      alpha=rfrf/rr;
      
      //     calculate p'=r'+alpha*p
      {
	for(int i=0;i<loc_vol;i++)
          {
	    bgp_load_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,r[i]);
	    bgp_load_spincolor(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,p[i]);
	    bgp_summassign_spincolor_prod_real(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,alpha);
	    bgp_save_spincolor(p[i],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
          }
      }

      //     calculate alphas=alpha*zfs*betas/zas*beta
      for(int imass=0;imass<nmass;imass++) alphas[imass]=alpha*zfs[imass]*betas[imass]/(zas[imass]*betaa);
      //     calculate ps'=zfs*r'+alphas*ps
      {
	for(int i=0;i<loc_vol;i++)
	  {
	    bgp_load_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,r[i]);
	    for(int imass=0;imass<nmass;imass++) 
	      if(flag[imass]==1)
		{
		  bgp_load_spincolor(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,ps[imass][i]);
		  bgp_assign_spincolor_prod_real(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,alphas[imass]);
		  bgp_summassign_spincolor_prod_real(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,zfs[imass]);
		  bgp_save_spincolor(ps[imass][i],B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32);
		}
	  }
      }
	    
	//     check over residual
	//if(rr*zfs[imass]<residue) flag[imass]=0;
        
	//     passes all f into actuals
	for(int imass=0;imass<nmass;imass++)
	  {
	    zps[imass]=zas[imass];
	    zas[imass]=zfs[imass];
	  }
      rr=rfrf;
      iter++;

      if(rank==0 && debug) printf("cgmms iter %d residue %g\n",iter,rfrf);
    }
  while(rfrf>stopping_residue && iter<niter);

  for(int imass=0;imass<nmass;imass++)  communicate_lx_spincolor_borders(sol[imass]);
  
  for(int imass=0;imass<nmass;imass++) free(ps[imass]);
  free(ps);
  free(s);
  free(p);
  free(r);
  free(t);
}
