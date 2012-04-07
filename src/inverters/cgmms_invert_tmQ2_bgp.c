#pragma once

void inv_tmQ2_cgmms_RL(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter,double st_res,double st_minres,int st_crit,int RL,spincolor *source)
{
  bgp_complex A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32;
  bgp_complex B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32;
  bgp_complex N0,N1,N2;
  
  double zps[nmass],zas[nmass],zfs[nmass],betas[nmass],alphas[nmass];
  double rr,rfrf,pap,alpha;
  double betap,betaa;
  int iter;
  int run_flag[nmass],nrun_mass=nmass;
  double final_res[nmass];
  double source_norm;
  
  spincolor *t=nissa_malloc("DD_temp",loc_vol+loc_bord,spincolor);
  spincolor *s=nissa_malloc("s",loc_vol,spincolor);
  spincolor *r=nissa_malloc("r",loc_vol,spincolor);
  spincolor *p=nissa_malloc("p",loc_vol+loc_bord,spincolor);
  spincolor *ps[nmass];
  for(int imass=0;imass<nmass;imass++) ps[imass]=nissa_malloc("ps",loc_vol,spincolor);

  //sol[*]=0
  for(int imass=0;imass<nmass;imass++)
    {
      memset(sol[imass],0,sizeof(spincolor)*loc_vol);
      run_flag[imass]=1;
      set_borders_invalid(sol[imass]);
    }
  
  //     -p=source
  //     -r=source
  //     -calculate Rr=(r,r)
  {
    complex cloc_rr={0,0};
    bgp_color_put_to_zero(N0,N1,N2);
    
    nissa_loc_vol_loop(ivol)
      {
	bgp_spincolor_load(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,source[ivol]);
	bgp_spincolor_save(p[ivol],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	bgp_spincolor_save(r[ivol],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	bgp_summassign_color_square_spincolor(N0,N1,N2,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
      }
    bgp_square_norm_color(N0,N1,N2);
    bgp_complex_save(cloc_rr,N0);
    set_borders_invalid(p);
    
    MPI_Allreduce(cloc_rr,&rr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    source_norm=rr;
    master_printf("Source norm: %lg\n",source_norm);
    
    master_printf("cgmms iter 0 rel. residues: ");
    for(int imass=0;imass<nmass;imass++) master_printf("%1.4e  ",1.0);
    master_printf("\n");
  }
  
  //     -betaa=1
  betaa=1;
  
  //     -ps=source
  {
    nissa_loc_vol_loop(ivol)
      {
	bgp_spincolor_load(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,source[ivol]);
	for(int imass=0;imass<nmass;imass++)
	  bgp_spincolor_save(ps[imass][ivol],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
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
      apply_tmQ2_RL(s,conf,kappa,m[0],t,RL,p);
      //     -pap=(p,s)=(p,Ap)
      {
	complex cloc_pap={0,0};
	bgp_color_put_to_zero(N0,N1,N2);
	
	nissa_loc_vol_loop(ivol)
	  {
	    bgp_spincolor_load(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,s[ivol]);
	    bgp_spincolor_load(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,p[ivol]);
	    bgp_summassign_color_realscalarprod_spincolor(N0,N1,N2,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32);
	  }
	bgp_square_norm_color(N0,N1,N2);
	bgp_complex_save(cloc_pap,N0);
	
	MPI_Allreduce(cloc_pap,&pap,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
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
		nissa_loc_vol_loop(ivol)
		  {
		    bgp_spincolor_load(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,sol[imass][ivol]);
		    bgp_spincolor_load(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,ps[imass][ivol]);
		    bgp_subtassign_spincolor_prod_double(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,betas[imass]);
		    bgp_spincolor_save(sol[imass][ivol],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
		  }
	      }
	      set_borders_invalid(sol[imass]);
	    }
	}
      
      //     calculate
      //     -r'=r+betaa*s=r+beta*Ap
      //     -rfrf=(r',r')
      {
	complex cloc_rfrf={0,0};
	bgp_color_put_to_zero(N0,N1,N2);
	
	nissa_loc_vol_loop(ivol)
	  {
	    bgp_spincolor_load(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,r[ivol]);
	    bgp_spincolor_load(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,s[ivol]);
	    bgp_summassign_spincolor_prod_double(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,betaa);
	    bgp_summassign_color_square_spincolor(N0,N1,N2,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	    bgp_spincolor_save(r[ivol],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	  }
	bgp_square_norm_color(N0,N1,N2);
	bgp_complex_save(cloc_rfrf,N0);
	
	MPI_Allreduce(cloc_rfrf,&rfrf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
      
      //     calculate alpha=rfrf/rr=(r',r')/(r,r)
      alpha=rfrf/rr;
      
      //     calculate p'=r'+alpha*p
      {
	nissa_loc_vol_loop(ivol)
	  {
	    bgp_spincolor_load(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,r[ivol]);
	    bgp_spincolor_load(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,p[ivol]);
	    bgp_summassign_spincolor_prod_double(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,alpha);
	    bgp_spincolor_save(p[ivol],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	  }
	set_borders_invalid(p);
      }
      
      for(int imass=0;imass<nmass;imass++)
	if(run_flag[imass]==1) //     calculate alphas=alpha*zfs*betas/zas*beta
	  alphas[imass]=alpha*zfs[imass]*betas[imass]/(zas[imass]*betaa);
      
      nissa_loc_vol_loop(ivol)
	{
	  bgp_spincolor_load(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,r[ivol]);
	  for(int imass=0;imass<nmass;imass++)
	    if(run_flag[imass]==1) //     calculate ps'=zfs*r'+alphas*ps
	      {
		bgp_spincolor_load(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,ps[imass][ivol]);
		bgp_assign_spincolor_prod_double(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,alphas[imass]);
		bgp_summassign_spincolor_prod_double(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,zfs[imass]);
		bgp_spincolor_save(ps[imass][ivol],B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32);
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
      nrun_mass=check_cgmms_residue_tmQ2_RL(run_flag,final_res,nrun_mass,rr,zfs,st_crit,st_res,st_minres,iter,nmass,m,source,conf,kappa,s,t,source_norm,RL,sol);
    }
  while(nrun_mass>0 && iter<niter);
  
  //print the final true residue
  
  for(int imass=0;imass<nmass;imass++)
    {
      double res,w_res,weight,max_res;
      apply_tmQ2_RL(s,conf,kappa,m[imass],t,RL,sol[imass]);
      {
	complex cloc_res={0,0};
	double locw_res=0;
	double locmax_res=0;
	double loc_weight=0;
	
	bgp_color_put_to_zero(N0,N1,N2);
	
	nissa_loc_vol_loop(ivol)
	  {
	    bgp_spincolor_load(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,s[ivol]);
	    bgp_spincolor_load(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,source[ivol]);
	    bgp_subtassign_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32,B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32);
	    bgp_squareassign_spincolor(A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	    bgp_summassign_color_spincolor(N0,N1,N2,A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	    
	    //calculate the weighted res
	    bgp_spincolor_save(s[ivol],A00,A01,A02,A10,A11,A12,A20,A21,A22,A30,A31,A32);
	    bgp_spincolor_load(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32,sol[imass][ivol]);
	    bgp_squareassign_spincolor(B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32);
	    bgp_spincolor_save(p[ivol],B00,B01,B02,B10,B11,B12,B20,B21,B22,B30,B31,B32);
	    set_borders_invalid(p);
	    
	    for(int id=0;id<4;id++)
	      for(int ic=0;ic<3;ic++)
		{
		  double plain_res=s[ivol][id][ic][0]+s[ivol][id][ic][1];
		  double point_weight=1/(p[ivol][id][ic][0]+p[ivol][id][ic][1]);
		  
		  locw_res+=plain_res*point_weight;
		  loc_weight+=point_weight;
		  if(plain_res>locmax_res) locmax_res=plain_res;
		}
	  }
	bgp_square_norm_color(N0,N1,N2);
	bgp_complex_save(cloc_res,N0);
	
	MPI_Reduce(cloc_res,&res,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&locw_res,&w_res,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&loc_weight,&weight,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&locmax_res,&max_res,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	
	w_res=w_res/weight*12*glb_vol;
	max_res*=12*glb_vol;
	
	master_printf("imass %d, rel residue true=%g approx=%g weighted=%g max=%g\n",imass,res/source_norm,final_res[imass]/source_norm,w_res,max_res);
      }
    }  
  
  for(int imass=0;imass<nmass;imass++) nissa_free(ps[imass]);
  nissa_free(s);
  nissa_free(p);
  nissa_free(r);
  nissa_free(t);  
}

