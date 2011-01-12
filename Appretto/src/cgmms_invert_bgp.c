#pragma once

#include "dirac_operator.c"
#include "su3.c"

double calculate_weighted_residue(spincolor *source,spincolor *sol,quad_su3 *conf,double kappa,double m,spincolor *s,spincolor *t,int dinf)
{
  apply_Q2(s,sol,conf,kappa,m,t,NULL,NULL);
  
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

int check_cgmms_residue(int *run_flag,double *residue_mass,int nrun,double rr,double *zfs,int st_crit,double st_res,double st_res2,int iter,spincolor **sol,int nmass,double *m,spincolor *source,quad_su3 *conf,double kappa,spincolor *s,spincolor *t)
{
  const int each=10;

  for(int imass=0;imass<nmass;imass++)
    if(run_flag[imass])
      {
	int fini=0;

	if(st_crit==sc_standard||st_crit==sc_unilevel)
	  {
	    residue_mass[imass]=rr*zfs[imass]*zfs[imass];
	    
	    if(st_crit==sc_standard) fini=(residue_mass[imass]<st_res2 || residue_mass[0]<st_res);
	    else fini=residue_mass[imass]<st_res;
	  }
	else if(st_crit==sc_weighted_norm2||st_crit==sc_weighted_norm_inf)
	  if(iter%each==0)
	    {  //locally weighted norm
	      communicate_lx_spincolor_borders(sol[imass]);
	      if(st_crit==sc_weighted_norm2)
		residue_mass[imass]=calculate_weighted_residue(source,sol[imass],conf,kappa,m[imass],s,t,2);
	      else residue_mass[imass]=calculate_weighted_residue(source,sol[imass],conf,kappa,m[imass],s,t,-1);
	
	      if(residue_mass[imass]<st_res) fini=1;
	    }
	
	if(fini)
	  {
	    run_flag[imass]=0;
	    nrun--;
	  }
      } 
  
  if(rank==0 && iter%each==0 && debug)
    {
      printf("cgmms iter %d residues: ",iter);
      for(int imass=0;imass<nmass;imass++) if(run_flag[imass]) printf("%1.4e  ",residue_mass[imass]);
      printf("\n");
    }
  
  return nrun;
}

void inv_Q2_cgmms(spincolor **sol,spincolor *source,spincolor **guess,quad_su3 *conf,double kappa,double *m,int nmass,int niter,double st_res,int st_crit)
{
  static double _Complex A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32;
  static double _Complex B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32;
  static double _Complex N0,N1,N2;
  static double _Complex R;

  double zps[nmass],zas[nmass],zfs[nmass],betas[nmass],alphas[nmass];
  double rr,rfrf,pap,alpha;
  double betap,betaa;
  int iter;
  int run_flag[nmass],nrun_mass=nmass;
  double final_res[nmass];
  double st_res2;

  spincolor *t=allocate_spincolor(loc_vol+loc_bord,"temporary for internal calculation of DD");
  spincolor *s=allocate_spincolor(loc_vol,"s in cgmms");
  spincolor *r=allocate_spincolor(loc_vol,"r in cgmms");
  spincolor *p=allocate_spincolor(loc_vol+loc_bord,"p in cgmms");
  spincolor **ps=(spincolor**)malloc(sizeof(spincolor*)*nmass);
  for(int imass=0;imass<nmass;imass++) ps[imass]=allocate_spincolor(loc_vol,"ps in cgmms");

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

    if(st_crit==sc_standard||st_crit==sc_unilevel) st_res*=rr;
    if(st_crit==sc_standard) st_res2=st_res*st_res;
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
      apply_Q2(s,p,conf,kappa,m[0],t,NULL,NULL);
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
          if(run_flag[imass]==1)
            {
	      zfs[imass]=zas[imass]*betap/(betaa*alpha*(1-zas[imass]/zps[imass])+betap*(1-(m[imass]-m[0])*(m[imass]+m[0])*betaa));
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

      for(int imass=0;imass<nmass;imass++)
        if(run_flag[imass]==1) //     calculate alphas=alpha*zfs*betas/zas*beta
	  alphas[imass]=alpha*zfs[imass]*betas[imass]/(zas[imass]*betaa);

      for(int i=0;i<loc_vol;i++)
	{
	  bgp_load_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,r[i]);
	  for(int imass=0;imass<nmass;imass++)
	    if(run_flag[imass]==1) //     calculate ps'=zfs*r'+alphas*ps
	      {
		bgp_load_spincolor(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,ps[imass][i]);
		bgp_assign_spincolor_prod_real(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,alphas[imass]);
		bgp_summassign_spincolor_prod_real(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,zfs[imass]);
		bgp_save_spincolor(ps[imass][i],B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32);
	      }
	}
  
      for(int imass=0;imass<nmass;imass++)
	if(run_flag[imass]==1) //     passes all f into actuals
	  {
	    zps[imass]=zas[imass];
	    zas[imass]=zfs[imass];
	  }

	rr=rfrf;
	iter++;
	
	//     check over residual
	nrun_mass=check_cgmms_residue(run_flag,final_res,nrun_mass,rr,zfs,st_crit,st_res,st_res2,iter,sol,nmass,m,source,conf,kappa,s,t);
    }
  while(nrun_mass>0 && iter<niter);
  
  for(int imass=0;imass<nmass;imass++)  communicate_lx_spincolor_borders(sol[imass]);
  
  //print the final true residue
  
  for(int imass=0;imass<nmass;imass++)
    {
      double res,w_res,weight,max_res;
      apply_Q2(s,sol[imass],conf,kappa,m[imass],t,NULL,NULL);
      {
	complex cloc_res={0,0};
	double locw_res=0;
	double locmax_res=0;
	double loc_weight=0;

        bgp_color_put_to_zero(N0,N1,N2);

	for(int i=0;i<loc_vol;i++)
	  {
	    bgp_load_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,s[i]);
            bgp_load_spincolor(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,source[i]);
	    bgp_subtassign_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32);
	    bgp_squareassign_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
            bgp_summassign_color_spincolor(N0,N1,N2,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);

	    //calculate the weighted res
	    bgp_save_spincolor(s[i],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
            bgp_load_spincolor(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,sol[imass][i]);
	    bgp_squareassign_spincolor(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32);
	    bgp_save_spincolor(p[i],B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32);
	    
	    for(int id=0;id<4;id++)
	      for(int ic=0;ic<3;ic++)
		{
		  double plain_res=s[i][id][ic][0]+s[i][id][ic][1];
		  double point_weight=1/(p[i][id][ic][0]+p[i][id][ic][1]);

		  locw_res+=plain_res*point_weight;
		  loc_weight+=point_weight;
		  if(plain_res>locmax_res) locmax_res=plain_res;
		}
	  }
        bgp_square_norm_color(N0,N1,N2);
	bgp_save_complex(cloc_res,N0);

	if(rank_tot>0)
	  {
	    MPI_Reduce(cloc_res,&res,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	    MPI_Reduce(&locw_res,&w_res,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	    MPI_Reduce(&loc_weight,&weight,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	    MPI_Reduce(&locmax_res,&max_res,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	  }
	else
	  {
	    res=cloc_res[0];
	    w_res=locw_res;
	    weight=loc_weight;
	    max_res=locmax_res;
	  }
	w_res=w_res/weight*12*glb_vol;
	max_res*=12*glb_vol;
	
	if(rank==0) printf("imass %d, residue true=%g approx=%g weighted=%g max=%g\n",imass,res,final_res[imass],w_res,max_res);
      }
    }  
  
  for(int imass=0;imass<nmass;imass++) free(ps[imass]);
  free(ps);
  free(s);
  free(p);
  free(r);
  free(t);
}
