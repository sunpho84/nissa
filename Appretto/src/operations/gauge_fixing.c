#pragma once

//determine the gauge transformation bringing to temporal gauge with T-1 timeslice diferent from id
void find_temporal_gauge_fixing_matr(su3 *fixm,quad_su3 *u)
{
  int loc_slice_area=loc_size[1]*loc_size[2]*loc_size[3];
  su3 *buf=NULL;
  
  //if the number of processors in the 0 dir is greater than 1 allocate room for border
  if(nproc_dir[0]>1) buf=appretto_malloc("buf",loc_slice_area,su3);

  //if we are on first proc slice put to identity the t=0 slice, otherwise receive it from previous proc slice
  if(proc_coord[0]==0)
    {
      for(int ivol=0;ivol<loc_vol;ivol++)
	if(glb_coord_of_loclx[ivol][0]==0)
	  su3_put_to_id(fixm[ivol]);
    }
  else
    if(nproc_dir[0]>1)
      MPI_Recv((void*)fixm,loc_slice_area,MPI_SU3,rank_neighdw[0],252,cart_comm,MPI_STATUS_IGNORE);
  
  //now go ahead along t
  int c[4];
  //loop over spatial slice
  for(c[1]=0;c[1]<loc_size[1];c[1]++)
    for(c[2]=0;c[2]<loc_size[2];c[2]++)
      for(c[3]=0;c[3]<loc_size[3];c[3]++)
	{
	  //bulk
	  for(c[0]=1;c[0]<loc_size[0];c[0]++)
	    {
	      int icurr=loclx_of_coord(c);
	      c[0]--;int iback=loclx_of_coord(c);c[0]++;
	      
	      su3_prod_su3(fixm[icurr],fixm[iback],u[iback][0]);
	    }
	  //border
	  if(nproc_dir[0]>1)
	    {
	      c[0]=loc_size[0]-1;int iback=loclx_of_coord(c);
	      c[0]=0;int icurr=loclx_of_coord(c);
	      
	      su3_prod_su3(buf[icurr],fixm[iback],u[iback][0]);
	    }
	  
	}
  
  //if we are not on last slice of processor send g to next slice
  if(proc_coord[0]!=(nproc_dir[0]-1) && nproc_dir[0]>1)
    MPI_Send((void*)buf,loc_slice_area,MPI_SU3,rank_neighup[0],252,cart_comm);

  if(nproc_dir[0]>1) appretto_free(buf);
}

//apply a gauge transformation to the conf
void gauge_transform_conf(quad_su3 *uout,su3 *g,quad_su3 *uin)
{
  if(rank_tot>1)
    {
      communicate_su3_borders(g);
      communicate_gauge_borders(uin);
    }
  
  su3 temp;
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int mu=0;mu<4;mu++)
      {
	su3_prod_su3_dag(temp,uin[ivol][mu],g[loclx_neighup[ivol][mu]]);
	su3_prod_su3(uout[ivol][mu],g[ivol],temp);
      }
}

//find the number of dir relevantfor the landau or coulomb condition
int find_nmu_gauge_cond(enum gauge_cond_type GT)
{
  int nmu=0;
  
  switch(GT)
    {
    case LANDAU: nmu=4;break;
    case COULOMB: nmu=3;break;
    default:crash("Calling compute_landau_delta with wrong gauge condit. %d (possible: %d or %d)",GT,LANDAU,COULOMB);break;
    }
  
  return nmu;
}

//overrelax the transformation
void overrelax(su3 *out,su3 *in,int N,double omega)
{
  //find coefficients
  double coef[N+1];
  for(int n=0;n<=N;n++) coef[n]=exp(lgamma(omega+1)-lgamma(omega+1-n)-lfact(n));
  
  //loop over volume
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      su3 t,o,w;

      su3_summ_real(w,in[ivol],-1); //subtract 1 from in[ivol]
      su3_copy(t,w);                //init t
      su3_put_to_id(o);             //output init
      
      //compute various powers
      for(int n=1;n<=N;n++)
	{ //summ 
	  for(int i=0;i<18;i++) ((double*)o)[i]+=coef[n]*((double*)t)[i];
	  safe_su3_prod_su3(t,t,w); //next power of w
	}
      
      //unitarize
      su3_unitarize(out[ivol],o);
    }
}

//apply the fast fourier acceleration
void fast_fourier_accelerate_fixing(su3 *g)
{  
  fft4d((complex*)g,(complex*)g,9,1);
  
  //apply the conditioning
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      //compute p2
      double p2=0;
      for(int mu=0;mu<4;mu++)
	{
	  int ix=glb_coord_of_loclx[ivol][mu];
	  double pmu=M_PI*ix/glb_size[mu];
	  double sinpmu=sin(pmu);
	  p2+=sinpmu*sinpmu;
	}
      p2/=4;
      
      //apply
      if(p2!=0)
	for(int i=0;i<3;i++)
	  for(int j=0;j<3;j++)
	    for(int ri=0;ri<2;ri++)
	      g[ivol][i][j][ri]/=p2;
    }
  
    fft4d((complex*)g,(complex*)g,9,-1);
}

//compute delta to fix landau or coulomb gauges
void compute_fixing_delta(su3 g,quad_su3 *conf,int ivol,enum gauge_cond_type GT)
{
  //reset g
  memset(g,0,sizeof(su3));
  
  //decide if landau or coulomb
  int nmu=find_nmu_gauge_cond(GT);
  
  //Calculate \sum_mu(U_mu(x-mu)-U_mu(x-mu)^dag-U_mu(x)+U^dag_mu(x))
  //and subtract the trace. It is computed just as the TA (traceless anti-herm)
  //this has 8 independent element (A+F+I)=0
  // ( 0,A) ( B,C) (D,E)
  // (-B,C) ( 0,F) (G,H)
  // (-D,E) (-G,H) (0,I)
  for(int mu=0;mu<nmu;mu++)
    {
      int b=loclx_neighdw[ivol][mu];
      
      //off-diagonal real parts
      g[0][1][0]+=conf[b][mu][0][1][0]-conf[b][mu][1][0][0]-conf[ivol][mu][0][1][0]+conf[ivol][mu][1][0][0]; //B
      g[0][2][0]+=conf[b][mu][0][2][0]-conf[b][mu][2][0][0]-conf[ivol][mu][0][2][0]+conf[ivol][mu][2][0][0]; //D
      g[1][2][0]+=conf[b][mu][1][2][0]-conf[b][mu][2][1][0]-conf[ivol][mu][1][2][0]+conf[ivol][mu][2][1][0]; //G
      
      //off diagonal imag parts
      g[0][1][1]+=conf[b][mu][0][1][1]+conf[b][mu][1][0][1]-conf[ivol][mu][0][1][1]-conf[ivol][mu][1][0][1]; //C
      g[0][2][1]+=conf[b][mu][0][2][1]+conf[b][mu][2][0][1]-conf[ivol][mu][0][2][1]-conf[ivol][mu][2][0][1]; //E
      g[1][2][1]+=conf[b][mu][1][2][1]+conf[b][mu][2][1][1]-conf[ivol][mu][1][2][1]-conf[ivol][mu][2][1][1]; //H
      
      //digonal imag parts
      g[0][0][1]+=conf[b][mu][0][0][1]+conf[b][mu][0][0][1]-conf[ivol][mu][0][0][1]-conf[ivol][mu][0][0][1]; //A
      g[1][1][1]+=conf[b][mu][1][1][1]+conf[b][mu][1][1][1]-conf[ivol][mu][1][1][1]-conf[ivol][mu][1][1][1]; //F
      g[2][2][1]+=conf[b][mu][2][2][1]+conf[b][mu][2][2][1]-conf[ivol][mu][2][2][1]-conf[ivol][mu][2][2][1]; //I
    }
  
  //compute the trace
  double T3=(g[0][0][1]+g[1][1][1]+g[2][2][1])/3;
  
  //subtract 1/3 of the trace from each element, so to make traceless the matrix
  g[0][0][1]-=T3;
  g[1][1][1]-=T3;
  g[2][2][1]-=T3;
  
  //fill the other parts
  
  //off-diagonal real parts
  g[1][0][0]=-g[0][1][0];
  g[2][0][0]=-g[0][2][0];
  g[2][1][0]=-g[1][2][0];
  
  //off diagonal imag parts
  g[1][0][1]=g[0][1][1];
  g[2][0][1]=g[0][2][1];
  g[2][1][1]=g[1][2][1];
}

//compute the steepest descent gauge fixing transformation
void find_steepest_descent_fixing(su3 *g,quad_su3 *conf,double alpha,enum gauge_cond_type GT)
{
  //loop over local sites
  for(int ivol=0;ivol<loc_vol;ivol++) compute_fixing_delta(g[ivol],conf,ivol,GT);
  
  fast_fourier_accelerate_fixing(g);

  //multiply by a/2 and add 1
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      for(int ic1=0;ic1<3;ic1++)
	{
	  for(int ic2=0;ic2<3;ic2++)
	    for(int ri=0;ri<2;ri++)
	      g[ivol][ic1][ic2][ri]*=alpha/2;
	  g[ivol][ic1][ic1][0]+=1;
	}
      su3_unitarize(g[ivol],g[ivol]);
    }
      
  //overrelax and unitarize
  //overrelax(g,g,3,1.75); //this is pointless
}

//here for future usage
double compute_gauge_fixing_functional(quad_su3 *conf,enum gauge_cond_type GT)
{  
  double loc_qual=0;
  double glb_qual;
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int mu=0;mu<4;mu++)
      for(int ic=0;ic<3;ic++)
	loc_qual+=conf[ivol][mu][ic][ic][0];
  
  //global reduction
  MPI_Allreduce(&loc_qual,&glb_qual,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return glb_qual/2/3/4/glb_vol;
}

//compute the quality of the gauge fixing
double compute_gauge_fixing_quality(quad_su3 *conf,enum gauge_cond_type GT)
{
  communicate_gauge_borders(conf);

  double loc_omega=0;
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      su3 delta;
      compute_fixing_delta(delta,conf,ivol,GT);
      
      //compute the trace of the square and summ it to omega
      for(int i=0;i<18;i++) loc_omega+=((double*)delta)[i]*((double*)delta)[i];
    }
  
  //global reduction
  double glb_omega;
  MPI_Allreduce(&loc_omega,&glb_omega,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return glb_omega/glb_vol/3;
}

//perform the landau or coulomb gauge fixing
void landau_or_coulomb_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double alpha,double precision,enum gauge_cond_type GT)
{
  memcpy(conf_out,conf_in,sizeof(quad_su3)*loc_vol);
  
  //fixing transformation
  su3 *fixm=appretto_malloc("fixm",loc_vol+loc_bord,su3);

  //fix iteratively up to reaching required precision
  int iter=0;
  double qual_out=compute_gauge_fixing_quality(conf_out,GT);
  master_printf("Iter 0 alpha %lg, quality: %lg (%lg req)\n",
		alpha,qual_out,precision);          
  do
    {
      //find the next fixing and compute its quality
      double tcomp=-take_time();
      find_steepest_descent_fixing(fixm,conf_out,alpha,GT);
      gauge_transform_conf(conf_out,fixm,conf_out);
      tcomp+=take_time();
      
      //compute change in quality
      double tqual=-take_time();
      double qual_in=qual_out;
      qual_out=compute_gauge_fixing_quality(conf_out,GT);
      double delta_qual=qual_out-qual_in;
      tqual+=take_time();
      
      iter++;
      master_printf("Iter %d alpha %lg, quality: %lg (%lg req) delta: %lg; %lg %lg s\n",
		    iter,alpha,qual_out,precision,delta_qual,tcomp,tqual);          
    }
  while(qual_out>=precision);
  
  master_printf("Final quality: %lg\n",qual_out);
  
  appretto_free(fixm);
}
//wrappers
void landau_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double alpha,double precision)
{landau_or_coulomb_gauge_fix(conf_out,conf_in,alpha,precision,LANDAU);}
void coulomb_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double alpha,double precision)
{landau_or_coulomb_gauge_fix(conf_out,conf_in,alpha,precision,COULOMB);}

