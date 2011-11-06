#pragma once

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

//overrelax the transformation
void overrelax(su3 out,su3 in,double w)
{
  double coef[5]={1,w,w*(w-1)/2,w*(w-1)*(w-2)/6,w*(w-1)*(w-2)*(w-3)/24};
  su3 f,t;
  
  su3_summ_real(f,in,-1);   //subtract 1 from in
  
  //ord 0
  su3_put_to_id(out);       //output init
  
  //ord 1
  su3_copy(t,f);
  su3_summ_the_prod_real(out,t,coef[1]);
  
  //ord 2-4
  for(int iord=2;iord<5;iord++)
    {
      safe_su3_prod_su3(t,t,f);
      su3_summ_the_prod_real(out,t,coef[iord]);
    }
}

//compute delta to fix landau or coulomb gauges
void compute_fixing_delta(su3 g,quad_su3 *conf,int ivol,enum gauge_cond_type GT)
{
}

//compute delta for the quality of landau or coulomb gauges
void compute_quality_delta(su3 g,quad_su3 *conf,int ivol,enum gauge_cond_type GT)
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

//horrible, horrifying routine of unknown meaning copied from APE
void exponentiate(su3 g,su3 a)
{
  //Exponentiate. Commented lines are original from APE (more or less)
  //up0=a[0][0]+a[1][1]~
  complex up0;
  complex_summ_conj2(up0,a[0][0],a[1][1]);
  //up1=a[0][1]-a[1][0]~
  complex up1;
  complex_subt_conj2(up1,a[0][1],a[1][0]);
  //icsi=1/sqrt(up0*up0~+up1*up1~)
  double icsi=1/sqrt(squared_complex_norm(up0)+squared_complex_norm(up1));
  //up0=up0*icsi
  complex_prod_real(up0,up0,icsi);
  //appuno=up0~
  complex appuno;
  complex_conj(appuno,up0);
  //up1=up1*icsi
  complex_prod_real(up1,up1,icsi);
  //appdue=-up1
  complex appdue;
  complex_prod_real(appdue,up1,-1);
  //n0=a[0][0]*appuno+a[1][0]*appdue
  complex n0;
  unsafe_complex_prod(n0,a[0][0],appuno);
  complex_summ_the_prod(n0,a[1][0],appdue);
  //n2=a[0][2]*appuno+a[1][2]*appdue
  complex n2;
  unsafe_complex_prod(n2,a[0][2],appuno);
  complex_summ_the_prod(n2,a[1][2],appdue);
  //n3=-a[0][0]*appdue~+a[1][0]*appuno~
  complex n3;
  unsafe_complex_conj2_prod(n3,a[1][0],appuno);
  complex_subt_the_conj2_prod(n3,a[0][0],appdue);
  //n6=a[2][0]
  complex n6={a[2][0][0],a[2][0][1]};
  //n5=-a[0][2]*appdue~+a[1][2]*appuno~
  complex n5;
  unsafe_complex_conj2_prod(n5,a[1][2],appuno);
  complex_subt_the_conj2_prod(n5,a[0][2],appdue);
  //n7=a[2][1]
  complex n7={a[2][1][0],a[2][1][1]};
  //up1=n5-n7~
  complex_subt_conj2(up1,n5,n7);
  //n4=-a[0][1]*appdue~+a[1][1]*appuno~
  complex n4;
  unsafe_complex_conj2_prod(n4,a[1][1],appuno);
  complex_subt_the_conj2_prod(n4,a[0][1],appdue);
  //n8=a[2][2]
  complex n8={a[2][2][0],a[2][2][1]};
  //up0=n8~+n4 
  complex_summ_conj1(up0,n8,n4);
  //icsi=1/sqrt(up0*up0~+up1*up1~)
  icsi=1/sqrt(squared_complex_norm(up0)+squared_complex_norm(up1));
  //up0=up0*icsi
  complex_prod_real(up0,up0,icsi);
  //bppuno=up0~
  complex bppuno;
  complex_conj(bppuno,up0);
  //up1=up1*icsi
  complex_prod_real(up1,up1,icsi);
  //bppdue=-up1
  complex bppdue={-up1[0],-up1[1]};
  //a[0][2]=n2
  a[0][2][0]=n2[0];
  a[0][2][1]=n2[1];
  //a[2][0]=-n3*bppdue~+n6*bppuno~
  unsafe_complex_conj2_prod(a[2][0],n6,bppuno);
  complex_subt_the_conj2_prod(a[2][0],n3,bppdue);
  //up1=a[0][2]-a[2][0]~
  complex_subt_conj2(up1,a[0][2],a[2][0]);
  //a[0][0]=n0
  a[0][0][0]=n0[0];
  a[0][0][1]=n0[1];
  //a[2][2]=-n5*bppdue~+n8*bppuno~
  unsafe_complex_conj2_prod(a[2][2],n8,bppuno);
  complex_subt_the_conj2_prod(a[2][2],n5,bppdue);
  //up0=a[2][2]~+a[0][0]
  complex_summ_conj1(up0,a[2][2],a[0][0]);
  //icsi=1/sqrt(up0*up0~+up1*up1~)
  icsi=1/sqrt(squared_complex_norm(up0)+squared_complex_norm(up1));
  //up0=up0*icsi
  complex_prod_real(up0,up0,icsi);
  //cppuno=up0~
  complex cppuno;
  complex_conj(cppuno,up0);
  //up1=up1*icsi
  complex_prod_real(up1,up1,icsi);
  //cppdue=-up1
  complex cppdue={-up1[0],-up1[1]};
  //e0=appuno
  complex e0={appuno[0],appuno[1]};
  //e1=appdue
  complex e1={appdue[0],appdue[1]};
  //e2=0
  //complex e2={0,0};
  //e3=-bppuno*appdue~
  complex e3={0,0};
  complex_subt_the_conj2_prod(e3,bppuno,appdue);
  //e4=bppuno*appuno~
  complex e4;
  unsafe_complex_conj2_prod(e4,bppuno,appuno);
  //e5=bppdue
  complex e5={bppdue[0],bppdue[1]};
  //e6=bppdue~*appdue~
  complex e6;
  unsafe_complex_conj_conj_prod(e6,bppdue,appdue);
  //e7=-bppdue~*appuno~
  complex e7={0,0};
  complex_subt_the_conj_conj_prod(e7,bppdue,appuno);
  //e8=bppuno~
  complex e8;
  complex_conj(e8,bppuno);
  //g[0][0]=cppuno*e0+cppdue*e6
  unsafe_complex_prod(g[0][0],cppuno,e0);
  complex_summ_the_prod(g[0][0],cppdue,e6);
  //g[0][1]=cppuno*e1+cppdue*e7
  unsafe_complex_prod(g[0][1],cppuno,e1);
  complex_summ_the_prod(g[0][1],cppdue,e7);
  //g[0][2]=cppdue*e8
  unsafe_complex_prod(g[0][2],cppdue,e8);
  //g[1][0]=e3
  g[1][0][0]=e3[0];
  g[1][0][1]=e3[1];
  //g[1][1]=e4
  g[1][1][0]=e4[0];
  g[1][1][1]=e4[1];
  //g[1][2]=e5
  g[1][2][0]=e5[0];
  g[1][2][1]=e5[1];
  //g[2][0]=-cppdue~*e0+cppuno~*e6
  unsafe_complex_conj1_prod(g[2][0],cppuno,e6);
  complex_subt_the_conj1_prod(g[2][0],cppdue,e0);
  //g[2][1]=-cppdue~*e1+cppuno~*e7
  unsafe_complex_conj1_prod(g[2][1],cppuno,e7);
  complex_subt_the_conj1_prod(g[2][1],cppdue,e1);
  //g[2][2]=cppuno~*e8
  unsafe_complex_conj1_prod(g[2][2],cppuno,e8);
}

//find the transformation bringing to the landau or coulomb gauge
void find_landau_or_coulomb_gauge_fixing_matr(su3 *fixm,quad_su3 *conf,double precision,enum gauge_cond_type GT)
{
  //allocate working conf
  quad_su3 *w_conf=appretto_malloc("Working conf",loc_vol+loc_bord,quad_su3);

  //set eo geometry, used to switch between different parity sites
  set_eo_geometry();
  
  //reset fixing transformation to unity
  for(int ivol=0;ivol<(loc_vol+loc_bord);ivol++) su3_put_to_id(fixm[ivol]);
  
  //fix iteratively up to reaching required precision
  int iter=0;
  double qual_out=compute_gauge_fixing_quality(conf_out,GT);
  master_printf("Iter 0, quality: %lg (%lg req)\n",
		qual_out,precision);
  
  //find the number of dir relevant for the landau or coulomb condition
  int nmu=0;
  switch(GT)
    {
    case LANDAU: nmu=4;break;
    case COULOMB: nmu=3;break;
    default:crash("Calling compute_landau_delta with wrong gauge condit. %d (possible: %d or %d)",GT,LANDAU,COULOMB);break;
    }
  
  //macro-loop in which the fixing is effectively applied to the original conf
  //this is done in order to avoid accumulated rounding errors
  do
    {
      //copy the gauge configuration on working fixing it with current transformation
      gauge_transform_conf(w_conf,fixm,conf);
      communicate_gauge_borders(w_conf);

      //loop fixing iteratively the working conf
      do
	{
	  //alternate even and odd
	  for(int par=0;par<2;par++)
	    {
	      // 1) first of all open
	      
	      //find the next fixing
	      
	      // 2) compute g=\sum_mu U_mu(x)+U^dag_mu(x-mu)

	      su3 g;
	      //first dir: reset and sum
	      int b=loclx_neighdw[ivol][0];
	      for(int ic1=0;ic1<3;ic1++)
		for(int ic2=0;ic2<3;ic2++)
		  {
		    g[ic1][ic2][0]=conf[ivol][0][ic1][ic2][0]+conf[b][0][ic2][ic1][0];
		    g[ic1][ic2][1]=conf[ivol][0][ic1][ic2][1]-conf[b][0][ic2][ic1][1];
		  }
	      //remaining dirs
	      for(int mu=1;mu<nmu;mu++)
		{
		  b=loclx_neighdw[ivol][mu];
		  
		  for(int ic1=0;ic1<3;ic1++)
		    for(int ic2=0;ic2<3;ic2++)
		      {
			g[ic1][ic2][0]+=conf[ivol][mu][ic1][ic2][0]+conf[b][mu][ic2][ic1][0];
			g[ic1][ic2][1]+=conf[ivol][mu][ic1][ic2][1]-conf[b][mu][ic2][ic1][1];
		      }
		}
	      
	      //exponentiate 
      double tcomp=-take_time();
      
      find_steepest_descent_fixing(fixm,conf_out,GT,iter);
      
      tcomp+=take_time();
      
      //compute change in quality
      double tqual=-take_time();
      double qual_in=qual_out;
      qual_out=compute_gauge_fixing_quality(conf_out,GT);
      double delta_qual=qual_out-qual_in;
      tqual+=take_time();
      
      iter++;
      master_printf("Iter %d, quality: %lg (%lg req) delta: %+3.3lg%%; %lg %lg s\n",
		    iter,qual_out,precision,delta_qual/(2*(qual_in+qual_out))*100,tcomp,tqual);          
    }
  while(qual_out>=precision);
  
  master_printf("Final quality: %lg\n",qual_out);
}

//compute the steepest descent gauge fixing transformation
void find_steepest_descent_fixing(su3 *g,quad_su3 *conf,enum gauge_cond_type GT,int iter)
{
  int par=iter%2;
  
  //loop over local sites
  for(int ivol=0;ivol<loc_vol;ivol++)
    if(loclx_parity[ivol]==par)
      {
	su3 a,b,c;
	compute_fixing_delta(a,conf,ivol,GT);
	exponentiate(b,a);
	overrelax(c,b,1.72);
	su3_unitarize(g[ivol],c);
      }
    else
      {
	memset(g[ivol],0,sizeof(su3));
	for(int ic=0;ic<3;ic++) g[ivol][ic][ic][0]=1;
      }
}

//compute the quality of the gauge fixing
double compute_gauge_fixing_quality(quad_su3 *conf,enum gauge_cond_type GT)
{
  communicate_gauge_borders(conf);

  double loc_omega=0;
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      su3 delta;
      compute_quality_delta(delta,conf,ivol,GT);
      
      //compute the trace of the square and summ it to omega
      for(int i=0;i<18;i++) loc_omega+=((double*)delta)[i]*((double*)delta)[i];
    }
  
  //global reduction
  double glb_omega;
  MPI_Allreduce(&loc_omega,&glb_omega,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return glb_omega/glb_vol/3;
}

//perform the landau or coulomb gauge fixing
void landau_or_coulomb_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double precision,enum gauge_cond_type GT)
{
  //allocate fixing matrix
  su3 *fixm=appretto_malloc("fixm",loc_vol+loc_bord,su3);
  
  //find fixing matrix
  find_landau_or_coulomb_gauge_fix(fixm,precision,GT);
  
  //apply the transformation
  gauge_transform_conf(conf_out,fixm,conf_out);  
  
  //free fixing matrix
  appretto_free(fixm);
}

//wrappers
void landau_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double precision)
{landau_or_coulomb_gauge_fix(conf_out,conf_in,precision,LANDAU);}
void coulomb_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double precision)
{landau_or_coulomb_gauge_fix(conf_out,conf_in,precision,COULOMB);}

