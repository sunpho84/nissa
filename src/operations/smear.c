#pragma once

//perform ape smearing
//be sure not to have border condition added
void ape_spatial_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep)
{
  quad_su3 *temp_conf=nissa_malloc("temp_conf",loc_vol+loc_bord+loc_edge,quad_su3);
  memcpy(smear_conf,origi_conf,sizeof(quad_su3)*loc_vol);
  
  if(debug_lvl) master_printf("APE smearing with alpha=%g, %d iterations\n",alpha,nstep);
      
  for(int istep=0;istep<nstep;istep++)
    {
      if(debug_lvl>1) master_printf("APE smearing with alpha=%g iteration %d of %d\n",alpha,istep,nstep);
      memcpy(temp_conf,smear_conf,sizeof(quad_su3)*loc_vol);
      
      //communicate the borders
      communicate_lx_quad_su3_borders(temp_conf);
      communicate_lx_gauge_edges(temp_conf);
      
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
		    unsafe_su3_prod_su3(temp1,temp_conf[A][nu],temp_conf[B][mu]);
		    unsafe_su3_prod_su3_dag(temp2,temp1,temp_conf[F][nu]);
		    su3_summ(stap,stap,temp2);
		        
		    int D=loclx_neighdw[A][nu];
		    int E=loclx_neighup[D][mu];
		    unsafe_su3_dag_prod_su3(temp1,temp_conf[D][nu],temp_conf[D][mu]);
		    unsafe_su3_prod_su3(temp2,temp1,temp_conf[E][nu]);
		    su3_summ(stap,stap,temp2);
		  }
	            
	      //create new link to be reunitarized
	      su3 prop_link;
	      for(int icol1=0;icol1<3;icol1++)
		for(int icol2=0;icol2<3;icol2++)
		  for(int ri=0;ri<2;ri++)
		    //prop_link[icol1][icol2][ri]=(1-alpha)*temp_conf[loc_site][mu][icol1][icol2][ri]+alpha/6*stap[icol1][icol2][ri];
		    prop_link[icol1][icol2][ri]=temp_conf[loc_site][mu][icol1][icol2][ri]+alpha*stap[icol1][icol2][ri];
	      
	      su3_unitarize_maximal_trace_projecting(smear_conf[loc_site][mu],prop_link);
	      //su3_unitarize_explicitly_inverting(smear_conf[loc_site][mu],prop_link);
	    }
	}
    }
      
  nissa_free(temp_conf);
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
  communicate_lx_spincolor_borders(smear_sc);
  memset(H,0,sizeof(spincolor)*loc_vol);
  
#ifndef BGP
  for(int l=0;l<loc_vol;l++)
    for(int id=0;id<4;id++)
      {
	for(int mu=1;mu<4;mu++)
	  {
	    int lup=loclx_neighup[l][mu];
	    int ldw=loclx_neighdw[l][mu];
	    
	    su3_summ_the_prod_color    (H[l][id],conf[l  ][mu],smear_sc[lup][id]);
	    su3_dag_summ_the_prod_color(H[l][id],conf[ldw][mu],smear_sc[ldw][id]);
	  }
	color_prod_double(H[l][id],H[l][id],kappa);
      }
#else
  bgp_complex H0,H1,H2;
  bgp_complex A0,A1,A2;
  bgp_complex B0,B1,B2;
  
  bgp_complex U00,U01,U02;
  bgp_complex U10,U11,U12;
  bgp_complex U20,U21,U22;

  bgp_complex V00,V01,V02;
  bgp_complex V10,V11,V12;
  bgp_complex V20,V21,V22;
  
  for(int l=0;l<loc_vol;l++)
    {
      for(int mu=1;mu<4;mu++)
	{
	  int lup=loclx_neighup[l][mu];
	  int ldw=loclx_neighdw[l][mu];
	  
	  bgp_cache_touch_su3(conf[l][mu]);
	  bgp_cache_touch_su3(conf[ldw][mu]);
	  bgp_cache_touch_spincolor(H[l]);
	  bgp_cache_touch_spincolor(smear_sc[lup]);
	  bgp_cache_touch_spincolor(smear_sc[ldw]);
	  
	  bgp_load_su3(U00,U01,U02,U10,U11,U12,U20,U21,U22,conf[l][mu]);
	  bgp_load_su3(V00,V01,V02,V10,V11,V12,V20,V21,V22,conf[ldw][mu]);
	  
	  for(int id=0;id<4;id++)
	    {
	      bgp_load_color(H0,H1,H2,H[l][id]);
	      bgp_load_color(A0,A1,A2,smear_sc[lup][id]);
	      bgp_load_color(B0,B1,B2,smear_sc[ldw][id]);
	      
	      bgp_summ_the_su3_prod_color(H0,H1,H2,U00,U01,U02,U10,U11,U12,U20,U21,U22,A0,A1,A2);
	      bgp_summ_the_su3_dag_prod_color(H0,H1,H2,V00,V01,V02,V10,V11,V12,V20,V21,V22,B0,B1,B2);

	      bgp_save_color(H[l][id],H0,H1,H2);
	    }
	}
      
      bgp_cache_touch_spincolor(H[l]);  
      for(int id=0;id<4;id++)
	{
	  bgp_load_color(A0,A1,A2,H[l][id]);  
	  bgp_color_prod_double(B0,B1,B2,A0,A1,A2,kappa);
	  bgp_save_color(H[l][id],B0,B1,B2);
	}
    }
  
#endif
}

//summ
void vol_spincolor_summassign(spincolor *smear_sc,spincolor *H)
{
  for(int l=0;l<loc_vol;l++)
    spincolor_summ(smear_sc[l],smear_sc[l],H[l]);
}

//prod with double
void vol_spincolor_prod_double(spincolor *out,spincolor *in,double r)
{
  for(int l=0;l<loc_vol;l++)
    spincolor_prod_double(out[l],in[l],r);
}

//jacobi smearing
void jacobi_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,double kappa,int niter)
{
  if(niter<1)
    {
      if(debug_lvl>1) master_printf("Skipping smearing (0 iter required)\n");
      if(smear_sc!=origi_sc) memcpy(smear_sc,origi_sc,sizeof(spincolor)*loc_vol);
    }
  else
    {
      spincolor *temp=nissa_malloc("temp",loc_vol+loc_bord,spincolor);//we do not know if smear_sc is allocated with bord
      spincolor *H=nissa_malloc("H",loc_vol+loc_bord,spincolor);
      double norm_fact=1/(1+6*kappa);
      communicate_lx_quad_su3_borders(conf);

      if(debug_lvl>1) master_printf("JACOBI smearing with kappa=%g, %d iterations\n",kappa,niter);
      
      //iter 0
      memcpy(temp,origi_sc,sizeof(spincolor)*loc_vol);
      
      //loop over jacobi iterations
      for(int iter=0;iter<niter;iter++)
	{
	  if(debug_lvl>1) master_printf("JACOBI smearing with kappa=%g iteration %d of %d\n",kappa,iter,niter);
	  
	  //apply kappa*H
	  smearing_apply_kappa_H(H,kappa,conf,temp);
#ifndef BGP
	  //add kappa*H
	  vol_spincolor_summassign(temp,H);
	  //dynamic normalization  
	  vol_spincolor_prod_double(temp,temp,norm_fact);
#else
	  bgp_complex A0,A1,A2;
	  bgp_complex B0,B1,B2;
	  bgp_complex C0,C1,C2;
	  
	  for(int l=0;l<loc_vol;l++)
	    {
	      bgp_cache_touch_spincolor(temp[l]);
	      bgp_cache_touch_spincolor(H[l]);
	      
	      for(int id=0;id<4;id++)
		{
		  bgp_load_color(A0,A1,A2,temp[l][id]);
		  bgp_load_color(B0,B1,B2,H[l][id]);
		  bgp_color_prod_double(C0,C1,C2,A0,A1,A2,norm_fact);
		  bgp_summassign_color_prod_double(C0,C1,C2,B0,B1,B2,norm_fact);
		  bgp_save_color(temp[l][id],C0,C1,C2);
		}
	    }
#endif
	}
      
      memcpy(smear_sc,temp,sizeof(spincolor)*loc_vol);
      
      nissa_free(H);
      nissa_free(temp);
    }
}

//smear a conf using hyp
//warning, the input conf needs to have edges allocate!
void hyp_smear_conf_dir(quad_su3 *sm_conf,quad_su3 *conf,double alpha0,double alpha1,double alpha2,int req_mu)
{
  //fill the dec2 remapping table
  int dec2_remap_index[4][4][4];
  int idec2_remap=0;
  for(int mu=0;mu<4;mu++)
    for(int nu=0;nu<4;nu++)
      for(int rho=0;rho<4;rho++)
	if(mu==nu || mu==rho || nu==rho || (req_mu>=0 && req_mu<=3 && mu!=req_mu && nu!=req_mu && rho!=req_mu)) dec2_remap_index[mu][nu][rho]=-1;
	else                                                                                                    dec2_remap_index[mu][nu][rho]=idec2_remap++;
  
  //fill the dec1 remapping table
  int dec1_remap_index[4][4];
  int idec1_remap=0;
  for(int mu=0;mu<4;mu++)
    for(int nu=0;nu<4;nu++)
      if(mu==nu || (req_mu>=0 && req_mu<=3 && mu!=req_mu && nu!=req_mu)) dec1_remap_index[mu][nu]=-1;
      else                                                               dec1_remap_index[mu][nu]=idec1_remap++;
  
  //communicate borders of original conf
  communicate_lx_quad_su3_borders(conf);
  communicate_lx_gauge_edges(conf);
  
  /////////////////////////////////////// second level decoration /////////////////////////////////
  
  //allocate dec2 conf
  su3 *dec2_conf[idec2_remap];
  for(int idec2=0;idec2<idec2_remap;idec2++) dec2_conf[idec2]=nissa_malloc("dec2_conf",loc_vol+loc_bord+loc_edge,su3);
  
  //loop over external index
  for(int mu=0;mu<4;mu++)
    //loop over the first decoration index
    for(int nu=0;nu<4;nu++)
      //loop over the second decoration index
      for(int rho=0;rho<4;rho++)
	if(nu!=mu && rho!=mu && rho!=nu && (req_mu<0 || req_mu>3 || mu==req_mu || nu==req_mu || rho==req_mu))
	  {
	    //find the remaining direction
	    int eta=0;
	    while(eta==mu || eta==nu || eta==rho) eta++;
	    
	    //find the remapped index
	    int ire0=dec2_remap_index[mu][nu][rho];
	    
	    //loop over local volume
	    for(int A=0;A<loc_vol;A++)
	      {
		//take original link
		su3 temp0;
		su3_prod_double(temp0,conf[A][mu],1-alpha2);
		
		//staple and temporary links
		su3 stap,temp1,temp2;
		
		//staple in the positive dir
		int B=loclx_neighup[A][eta];
		int F=loclx_neighup[A][mu];
		unsafe_su3_prod_su3(temp1,conf[A][eta],conf[B][mu]);
		unsafe_su3_prod_su3_dag(stap,temp1,conf[F][eta]);
		
		//staple in the negative dir
		int D=loclx_neighdw[A][eta];
		int E=loclx_neighup[D][mu];
		unsafe_su3_dag_prod_su3(temp1,conf[D][eta],conf[D][mu]);
		unsafe_su3_prod_su3(temp2,temp1,conf[E][eta]);
		su3_summ(stap,stap,temp2);
		
		//summ the two staples with appropriate coef
		su3_summ_the_prod_double(temp0,stap,alpha2/2);
		
		//project the resulting link onto su3
		su3_unitarize_maximal_trace_projecting(dec2_conf[ire0][A],temp0);
	      }
	    
	    //communicate borders for future usage
	    communicate_lx_su3_borders(dec2_conf[ire0]);
	    communicate_lx_su3_edges(dec2_conf[ire0]);
	  }
  
  /////////////////////////////////////// first level decoration /////////////////////////////////
  
  //allocate dec1 conf
  su3 *dec1_conf[idec1_remap];
  for(int idec1=0;idec1<idec1_remap;idec1++) dec1_conf[idec1]=nissa_malloc("dec1_conf",loc_vol+loc_bord+loc_edge,su3);
  
  //loop over external index
  for(int mu=0;mu<4;mu++)
    //loop over the first decoration index
    for(int nu=0;nu<4;nu++)
      if(nu!=mu && (req_mu<0 || req_mu>3 || mu==req_mu || nu==req_mu))
	//loop over local volume
	for(int A=0;A<loc_vol;A++)
	  {
	    //take original link
	    su3 temp0;
	    su3_prod_double(temp0,conf[A][mu],1-alpha1);
	    
	    //reset the staple
	    su3 stap;
	    memset(stap,0,sizeof(su3));
	    
	    //find the remapped index
	    int ire0=dec1_remap_index[mu][nu];
	    
	    //loop over the second decoration index
	    for(int rho=0;rho<4;rho++)
	      if(rho!=mu && rho!=nu)
		{
		  su3 temp1,temp2;
		  
		  //find the two remampped indices
		  int ire1=dec2_remap_index[rho][nu][mu];
		  int ire2=dec2_remap_index[mu][rho][nu];
		  
		  //staple in the positive dir
		  int B=loclx_neighup[A][rho];
		  int F=loclx_neighup[A][mu];
		  unsafe_su3_prod_su3(temp1,dec2_conf[ire1][A],dec2_conf[ire2][B]);
		  unsafe_su3_prod_su3_dag(temp2,temp1,dec2_conf[ire1][F]);
		  su3_summ(stap,stap,temp2);
		  
		  //staple in the negative dir
		  int D=loclx_neighdw[A][rho];
		  int E=loclx_neighup[D][mu];
		  unsafe_su3_dag_prod_su3(temp1,dec2_conf[ire1][D],dec2_conf[ire2][D]);
		  unsafe_su3_prod_su3(temp2,temp1,dec2_conf[ire1][E]);
		  su3_summ(stap,stap,temp2);
		  
		  //summ the two staples with appropriate coef
		  su3_summ_the_prod_double(temp0,stap,alpha1/4);
		  
		  //project the resulting link onto su3
		  su3_unitarize_maximal_trace_projecting(dec1_conf[ire0][A],temp0);
		}
	    
	    //communicate borders for future usage
	    communicate_lx_su3_borders(dec1_conf[ire0]);
	    communicate_lx_su3_edges(dec1_conf[ire0]);
	  }
  
  //free dec2
  for(int idec2=0;idec2<idec2_remap;idec2++) nissa_free(dec2_conf[idec2]);
  
  /////////////////////////////////////// zero level decoration /////////////////////////////////
  
  //loop over external index
  for(int mu=0;mu<4;mu++)
    //loop over local volume
    for(int A=0;A<loc_vol;A++)
      {
	//take original link
	su3 temp0;
	//exclude all dir apart from req_mu, if req_mu is in the range [0:3]
	if(req_mu>=0 && req_mu<=3 && mu!=req_mu)
	  {if(sm_conf!=conf) su3_copy(sm_conf[A][mu],conf[A][mu]);}
	else
	  {
	    su3_prod_double(temp0,conf[A][mu],1-alpha0);
	    
	    //reset the staple
	    su3 stap;
	    memset(stap,0,sizeof(su3));
	    
	    //loop over the first decoration index
	    for(int nu=0;nu<4;nu++)
	      if(nu!=mu)
		{
		  su3 temp1,temp2;
		  
		  //find the two remampped indices
		  int ire1=dec1_remap_index[nu][mu];
		  int ire2=dec1_remap_index[mu][nu];
		  
		  //staple in the positive dir
		  int B=loclx_neighup[A][nu];
		  int F=loclx_neighup[A][mu];
		  unsafe_su3_prod_su3(temp1,dec1_conf[ire1][A],dec1_conf[ire2][B]);
		  unsafe_su3_prod_su3_dag(temp2,temp1,dec1_conf[ire1][F]);
		  su3_summ(stap,stap,temp2);
		  
		  //staple in the negative dir
		  int D=loclx_neighdw[A][nu];
		  int E=loclx_neighup[D][mu];
		  unsafe_su3_dag_prod_su3(temp1,dec1_conf[ire1][D],dec1_conf[ire2][D]);
		  unsafe_su3_prod_su3(temp2,temp1,dec1_conf[ire1][E]);
		  su3_summ(stap,stap,temp2);
		  
		  //summ the two staples with appropriate coef
		  su3_summ_the_prod_double(temp0,stap,alpha0/6);
		  
		  //project the resulting link onto su3
		  su3_unitarize_maximal_trace_projecting(sm_conf[A][mu],temp0);
		}
	  }
      }
  
  //free dec1
  for(int idec1=0;idec1<idec1_remap;idec1++) nissa_free(dec1_conf[idec1]);
}

//hyp smear all the dirs
void hyp_smear_conf(quad_su3 *sm_conf,quad_su3 *conf,double alpha0,double alpha1,double alpha2)
{
  hyp_smear_conf_dir(sm_conf,conf,alpha0,alpha1,alpha2,-1);
}
