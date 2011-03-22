#pragma once

//perform ape smearing
//be sure not to have border condition added
void ape_smearing(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep)
{
  quad_su3 *temp_conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"temp_conf");
  memcpy(smear_conf,origi_conf,sizeof(quad_su3)*loc_vol);
  
  for(int istep=0;istep<nstep;istep++)
    {
      if(rank==0 && debug) printf("APE smearing with alpha=%g iteration %d of %d\n",alpha,istep,nstep);
      memcpy(temp_conf,smear_conf,sizeof(quad_su3)*loc_vol);
      
      //communicate the borders
      communicate_gauge_borders(temp_conf);
      communicate_gauge_edges(temp_conf);
      
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
	{     
	  for(int mu=1;mu<4;mu++)
	    {
	      //calculate staples
	      su3 stap,temp1,temp2;
	      memset(stap,0,sizeof(su3));
	      for(int nu=1;nu<4;nu++)                   //  E---F---C   
		if(nu!=mu)                              //  |   |   | mu
		  {                                     //  D---A---B   
		    int A=loc_site;                     //        nu    
		    int B=loclx_neighup[A][nu];
		    int F=loclx_neighup[A][mu];
		    su3_prod_su3(temp1,temp_conf[A][nu],temp_conf[B][mu]);
		    su3_prod_su3_dag(temp2,temp1,temp_conf[F][nu]);
		    su3_summ(stap,stap,temp2);
		        
		    int D=loclx_neighdw[A][nu];
		    int E=loclx_neighup[D][mu];
		    su3_dag_prod_su3(temp1,temp_conf[D][nu],temp_conf[D][mu]);
		    su3_prod_su3(temp2,temp1,temp_conf[E][nu]);
		    su3_summ(stap,stap,temp2);
		  }
	            
	      //create new link to be reunitarized
	      su3 prop_link;
	      for(int icol1=0;icol1<3;icol1++)
		for(int icol2=0;icol2<3;icol2++)
		  for(int ri=0;ri<2;ri++)
		    //prop_link[icol1][icol2][ri]=(1-alpha)*temp_conf[loc_site][mu][icol1][icol2][ri]+alpha/6*stap[icol1][icol2][ri];
		    prop_link[icol1][icol2][ri]=temp_conf[loc_site][mu][icol1][icol2][ri]+alpha*stap[icol1][icol2][ri];
	            
	      su3_unitarize(smear_conf[loc_site][mu],prop_link);
	    }
	}
    }
      
  free(temp_conf);
}

//return the spatial density profile at a fixed timeslice, at radius r from or_pos
void density_profile(double *glb_rho,spincolor *sp,int *or_pos)
{
  int L=glb_size[1];
  
  int glb_n[L],loc_n[L]; 
  double loc_rho[L];
  
  memset(loc_rho,0,sizeof(double)*L);
  memset(loc_n,0,sizeof(int)*L);

  for(int l=0;l<loc_vol;l++)
    if(glb_coord_of_loclx[l][0]==or_pos[0]||or_pos[0]<0)
      {
        //calculate distance
        double r2=0;
        int d;
        for(int mu=1;mu<4;mu++)
          {
            int x=glb_coord_of_loclx[l][mu]-or_pos[mu];
            if(x>=L/2) x-=L;
            r2+=x*x;
          }
        d=(int)sqrt(r2);

        //calculate the number of points at distance d
        loc_n[d]++;
        
        //calculate norm of the source
        for(int id=0;id<4;id++)
          for(int ic=0;ic<3;ic++)
            for(int ri=0;ri<2;ri++)
              loc_rho[d]+=sp[l][id][ic][ri]*sp[l][id][ic][ri];
      }
  
  //final reduction
  MPI_Reduce(loc_rho,glb_rho,L,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(loc_n,glb_n,L,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  
  //normalize
  for(int d=0;d<L;d++) if(glb_n[d]) glb_rho[d]=sqrt(glb_rho[d]/glb_n[d]);
}

//apply kappa*H to a spincolor
void smearing_apply_kappa_H(spincolor *H,double kappa,quad_su3 *conf,spincolor *smear_sc,int timeslice)
{
  memset(H,0,sizeof(spincolor)*loc_vol);
  for(int l=0;l<loc_vol;l++)
    if(glb_coord_of_loclx[l][0]==timeslice||timeslice<0)
      for(int id=0;id<4;id++)
	{
	  for(int mu=1;mu<4;mu++)
	    {
	      int lup=loclx_neighup[l][mu];
	      int ldw=loclx_neighdw[l][mu];
	      
	      su3_summ_the_color_prod    (H[l][id],conf[l  ][mu],smear_sc[lup][id]);
	      su3_dag_summ_the_color_prod(H[l][id],conf[ldw][mu],smear_sc[ldw][id]);
	    }
	  color_prod_real(H[l][id],H[l][id],kappa);
	}
}

//summ
void vol_spincolor_summassign(spincolor *smear_sc,spincolor *H,int timeslice)
{
  for(int l=0;l<loc_vol;l++)
    if(glb_coord_of_loclx[l][0]==timeslice||timeslice==-1)
      spincolor_summ(smear_sc[l],smear_sc[l],H[l]);
}

//normalize
void vol_assign_spincolor_prod_real(spincolor *sc,double c,int timeslice)
{
  for(int l=0;l<loc_vol;l++)
    if(glb_coord_of_loclx[l][0]==timeslice||timeslice==-1)
      assign_spincolor_prod_real(sc[l],c);
}

//normalize
void vol_spincolor_normalize(spincolor *smear_sc,int timeslice)
{
  //calculate global norm
  double loc_norm=0;
  for(int l=0;l<loc_vol;l++)
    if(glb_coord_of_loclx[l][0]==timeslice||timeslice==-1)
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    loc_norm+=smear_sc[l][id][ic][ri]*smear_sc[l][id][ic][ri];
  double norm;
  MPI_Allreduce(&loc_norm,&norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  norm=sqrt(norm);
  
  //normalize
  for(int l=0;l<loc_vol;l++)      
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  smear_sc[l][id][ic][ri]/=norm;
}

//jacobi smearing on the 
void jacobi_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,double kappa,int niter,int timeslice)
{
  spincolor *H_old=allocate_spincolor(loc_vol+loc_bord,"H_old");
  spincolor *H_new=allocate_spincolor(loc_vol+loc_bord,"H_new");
  
  communicate_gauge_borders(conf);
  
  //iter 0
  memcpy(smear_sc,origi_sc,sizeof(spincolor)*loc_vol);
  memcpy(H_new,origi_sc,sizeof(spincolor)*loc_vol);
  
  for(int iter=0;iter<niter;iter++)
    {
      memcpy(H_old,H_new,sizeof(spincolor)*loc_vol);
      communicate_lx_spincolor_borders(H_old);
      
      //apply kappa*H to H_old
      smearing_apply_kappa_H(H_new,kappa,conf,H_old,timeslice);
      
      //add kappa*H
      vol_spincolor_summassign(smear_sc,H_new,timeslice);
    }
  
  //final normalization  
  vol_spincolor_normalize(smear_sc,timeslice);
  
  free(H_old);
  free(H_new);
}

double vol_spincolor_norm(spincolor *smear_sc,int timeslice)
{
  double n2=0;
  for(int l=0;l<loc_vol;l++)
    if(glb_coord_of_loclx[l][0]==timeslice||timeslice==-1)
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    n2+=pow(smear_sc[l][id][ic][ri],2);
  
  return n2;
}  

//jacobi smearing on the 
void dina_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,double kappa,int niter,int timeslice)
{
  spincolor *H=allocate_spincolor(loc_vol+loc_bord,"H");
 
  communicate_gauge_borders(conf);
  
  //iter 0
  memcpy(smear_sc,origi_sc,sizeof(spincolor)*loc_vol);
  double n2_ini=vol_spincolor_norm(smear_sc,timeslice);

  //loop over jacobi iterations
  for(int iter=0;iter<niter;iter++)
    {
      communicate_lx_spincolor_borders(smear_sc);
      
      //apply kappa*H
      smearing_apply_kappa_H(H,kappa,conf,smear_sc,timeslice);

      //add kappa*H
      vol_spincolor_summassign(smear_sc,H,timeslice);
      
      //dynamic normalization  
      vol_assign_spincolor_prod_real(smear_sc,1/(1+6*kappa),timeslice);
    }
  
  double n2_fin=vol_spincolor_norm(smear_sc,timeslice);
  //final normalization  
  vol_assign_spincolor_prod_real(smear_sc,sqrt(n2_ini/n2_fin),timeslice);
  
  free(H);
}

