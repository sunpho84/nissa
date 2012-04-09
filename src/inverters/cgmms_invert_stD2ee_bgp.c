#pragma once

void inv_stD2ee_cgmm2s(color **sol,quad_su3 **conf,double *m2,int nmass,int niter,double st_res,double st_minres,int st_crit,color *source)
{
  bgp_complex A00,A01,A02;
  bgp_complex B00,B01,B02;
  bgp_complex N0,N1,N2;
  
  double zps[nmass],zas[nmass],zfs[nmass],betas[nmass],alphas[nmass];
  double rr,rfrf,pap,alpha;
  double betap,betaa;
  int iter;
  int run_flag[nmass],nrun_mass=nmass;
  double final_res[nmass];
  double source_norm;
  
  color *t=nissa_malloc("DD_temp",loc_volh+loc_bordh,color);
  color *s=nissa_malloc("s",loc_volh,color);
  color *r=nissa_malloc("r",loc_volh,color);
  color *p=nissa_malloc("p",loc_volh+loc_bordh,color);
  color *ps[nmass];
  for(int imass=0;imass<nmass;imass++) ps[imass]=nissa_malloc("ps",loc_volh,color);
  
  //sol[*]=0
  for(int imass=0;imass<nmass;imass++)
    {
      memset(sol[imass],0,sizeof(color)*loc_volh);
      run_flag[imass]=1;
      set_borders_invalid(sol[imass]);
    }
  
  //     -p=source
  //     -r=source
  //     -calculate Rr=(r,r)
  {
    complex cloc_rr={0,0};
    bgp_color_put_to_zero(N0,N1,N2);
    
    nissa_loc_volh_loop(ivol)
    {
      bgp_color_load(A00,A01,A02, source[ivol]);
      bgp_color_save(p[ivol], A00,A01,A02);
      bgp_color_save(r[ivol], A00,A01,A02);
      bgp_summassign_square_color(N0,N1,N2, A00,A01,A02);
    }
    bgp_square_norm_color(N0,N1,N2);
    bgp_complex_save(cloc_rr,N0);
    set_borders_invalid(p);
    
    MPI_Allreduce(cloc_rr,&rr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    source_norm=rr;
    master_printf("Source norm: %lg\n",source_norm);
    
    master_printf(" cgmms iter 0 rel. residues: ");
    for(int imass=0;imass<nmass;imass++) master_printf("%1.4e  ",1.0);
    master_printf("\n");
  }
  
  //     -betaa=1
  betaa=1;
  
  //     -ps=source
  {
    nissa_loc_volh_loop(ivol)
      {
	bgp_color_load(A00,A01,A02, source[ivol]);
	for(int imass=0;imass<nmass;imass++)
	  bgp_color_save(ps[imass][ivol], A00,A01,A02);
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
  int nrequest=0;
  MPI_Request request[16];
  do
    {
      //     -s=Ap
      if(nrequest!=0) finish_communicating_ev_color_borders(&nrequest,request,p);      
      apply_st2Doe(t,conf,p);
      apply_st2Deo(s,conf,t);
      
      //     -pap=(p,s)=(p,Ap)
      {
        complex cloc_pap={0,0};
        bgp_color_put_to_zero(N0,N1,N2);
        
        nissa_loc_volh_loop(ivol)
	{
	  bgp_color_load(A00,A01,A02, s[ivol]);
	  bgp_assign_color_prod_double(A00,A01,A02, 0.25);
	  bgp_color_save(s[ivol], A00,A01,A02);
	  bgp_color_load(B00,B01,B02, p[ivol]);
	  bgp_summassign_scalarprod_color(N0,N1,N2, A00,A01,A02, B00,B01,B02);
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
              zfs[imass]=zas[imass]*betap/(betaa*alpha*(1-zas[imass]/zps[imass])+betap*(1-m2[imass]*betaa));
              betas[imass]=betaa*zfs[imass]/zas[imass];
              
              {
                nissa_loc_volh_loop(ivol)
		{
		  bgp_color_load(A00,A01,A02, sol[imass][ivol]);
		  bgp_color_load(B00,B01,B02, ps[imass][ivol]);
		  bgp_subtassign_color_prod_double(A00,A01,A02, B00,B01,B02, betas[imass]);
		  bgp_color_save(sol[imass][ivol], A00,A01,A02);
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
        
        nissa_loc_volh_loop(ivol)
	{
	  bgp_color_load(A00,A01,A02, r[ivol]);
	  bgp_color_load(B00,B01,B02, s[ivol]);
	  bgp_summassign_color_prod_double(A00,A01,A02, B00,B01,B02, betaa);
	  bgp_summassign_square_color(N0,N1,N2, A00,A01,A02);
	  bgp_color_save(r[ivol], A00,A01,A02);
	}
        bgp_square_norm_color(N0,N1,N2);
        bgp_complex_save(cloc_rfrf,N0);
        
        MPI_Allreduce(cloc_rfrf,&rfrf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }

      //     calculate alpha=rfrf/rr=(r',r')/(r,r)
      alpha=rfrf/rr;
      
      //     calculate p'=r'+alpha*p
      {
        nissa_loc_volh_loop(ivol)
	{
	  bgp_color_load(A00,A01,A02, r[ivol]);
	  bgp_color_load(B00,B01,B02, p[ivol]);
	  bgp_summassign_color_prod_double(A00,A01,A02, B00,B01,B02, alpha);
	  bgp_color_save(p[ivol], A00,A01,A02);
	}
        set_borders_invalid(p);
	start_communicating_ev_color_borders(&nrequest,request,p);
      }

      for(int imass=0;imass<nmass;imass++)
        if(run_flag[imass]==1) //     calculate alphas=alpha*zfs*betas/zas*beta
          alphas[imass]=alpha*zfs[imass]*betas[imass]/(zas[imass]*betaa);
      
      nissa_loc_volh_loop(ivol)
      {
	bgp_color_load(A00,A01,A02, r[ivol]);
	for(int imass=0;imass<nmass;imass++)
	  if(run_flag[imass]==1) //     calculate ps'=zfs*r'+alphas*ps
	    {
	      bgp_color_load(B00,B01,B02, ps[imass][ivol]);
	      bgp_assign_color_prod_double(B00,B01,B02, alphas[imass]);
	      bgp_summassign_color_prod_double(B00,B01,B02, A00,A01,A02, zfs[imass]);
	      bgp_color_save(ps[imass][ivol], B00,B01,B02);
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
      nrun_mass=check_cgmm2s_residue_stD2ee(run_flag,final_res,nrun_mass,rr,zfs,st_crit,st_res,st_minres,iter,nmass,m2,source,conf,s,t,source_norm,sol);
    }
  while(nrun_mass>0 && iter<niter);
  
  //print the final true residue
  
  for(int imass=0;imass<nmass;imass++)
    {
      double res,w_res,weight,max_res;
      apply_stD2ee(s,conf,t,sqrt(m2[imass]),sol[imass]);
      {
        complex cloc_res={0,0};
        double locw_res=0;
        double locmax_res=0;
        double loc_weight=0;
        
        bgp_color_put_to_zero(N0,N1,N2);
        
        nissa_loc_volh_loop(ivol)
	{
	  bgp_color_load(A00,A01,A02, s[ivol]);
	  bgp_color_load(B00,B01,B02, source[ivol]);
	  bgp_subtassign_color(A00,A01,A02, B00,B01,B02);
	  bgp_squareassign_color(A00,A01,A02);
	  bgp_summassign_color(N0,N1,N2, A00,A01,A02);
            
	  //calculate the weighted res
	  bgp_color_save(s[ivol], A00,A01,A02);
	  bgp_color_load(B00,B01,B02, sol[imass][ivol]);
	  bgp_squareassign_color(B00,B01,B02);
	  bgp_color_save(p[ivol], B00,B01,B02);
	  set_borders_invalid(p);
	  
	  for(int ic=0;ic<3;ic++)
	    {
	      double plain_res=s[ivol][ic][0]+s[ivol][ic][1];
	      double point_weight=1/(p[ivol][ic][0]+p[ivol][ic][1]);
              
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
        
        w_res=w_res/weight*3*glb_volh;
        max_res*=3*glb_volh;
        
        master_printf("imass %d, rel residue true=%g approx=%g weighted=%g max=%g\n",imass,res/source_norm,final_res[imass]/source_norm,w_res,max_res);
      }
    }
  
  master_printf(" Total cgmms iterations: %d\n",iter);
  
  for(int imass=0;imass<nmass;imass++) nissa_free(ps[imass]);
  nissa_free(s);
  nissa_free(p);
  nissa_free(r);
  nissa_free(t);
}
