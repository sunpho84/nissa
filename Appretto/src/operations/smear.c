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
            if(x<-L/2) x+=L;
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
  for(int d=0;d<L;d++) if(glb_n[d])
    glb_rho[d]=sqrt(glb_rho[d]/glb_n[d]);
}

//apply kappa*H to a spincolor
void smearing_apply_kappa_H(spincolor *H,double kappa,quad_su3 *conf,spincolor *smear_sc)
{
  memset(H,0,sizeof(spincolor)*loc_vol);
  
  communicate_lx_spincolor_borders(smear_sc);
  
  for(int l=0;l<loc_vol;l++)
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
void vol_spincolor_summassign(spincolor *smear_sc,spincolor *H)
{
  for(int l=0;l<loc_vol;l++)
    spincolor_summ(smear_sc[l],smear_sc[l],H[l]);
}

//prod with real
void vol_assign_spincolor_prod_real(spincolor *sc,double c)
{
  for(int l=0;l<loc_vol;l++)
    assign_spincolor_prod_real(sc[l],c);
}

//jacobi smearing
void jacobi_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,double kappa,int niter)
{
  spincolor *H=allocate_spincolor(loc_vol+loc_bord,"H");
  
  communicate_gauge_borders(conf);
  
  //iter 0
  if(smear_sc!=origi_sc) memcpy(smear_sc,origi_sc,sizeof(spincolor)*loc_vol);
  
  //loop over jacobi iterations
  for(int iter=0;iter<niter;iter++)
    {
      //apply kappa*H
      smearing_apply_kappa_H(H,kappa,conf,smear_sc);

      //add kappa*H
      vol_spincolor_summassign(smear_sc,H);
      
      //dynamic normalization  
      vol_assign_spincolor_prod_real(smear_sc,1/(1+6*kappa));
    }
  
  free(H);
}

