#pragma once

//Local trace of the product of gamma1 * spinspin1^dag * gamma2 * spinspin2
void site_trace_g_sdag_g_s(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2)
{
  spinspin t1,t2;
  
  spinspin_dirac_spinspindag_prod(t1,g1,s1);
  spinspin_dirac_spinspin_prod(t2,g2,s2);

  trace_prod_spinspins(c,t1,t2);
}

//the eight
void site_trace_g_sdag_g_s_g_sdag_g_s(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2, dirac_matr *g3,spinspin s3,dirac_matr *g4,spinspin s4)
{
  spinspin t1,t2,t12,t3,t4,t34;
  
  spinspin_dirac_spinspindag_prod(t1,g1,s1);
  spinspin_dirac_spinspin_prod(t2,g2,s2);
  spinspin_spinspin_prod(t12,t1,t2);

  spinspin_dirac_spinspindag_prod(t3,g3,s3);
  spinspin_dirac_spinspin_prod(t4,g4,s4);
  spinspin_spinspin_prod(t34,t3,t4);
  
  trace_prod_spinspins(c,t12,t34);
}

//12 index prop
void site_trace_g_ccss_dag_g_ccss(complex c,dirac_matr *g1,su3spinspin s1,dirac_matr *g2,su3spinspin s2)
{
  c[0]=c[1]=0;
  //Color loop
  for(int ic1=0;ic1<3;ic1++)
    for(int ic2=0;ic2<3;ic2++)
      {
	spinspin t1,t2;
	
	spinspin_dirac_spinspindag_prod(t1,g1,s1[ic2][ic1]);
	spinspin_dirac_spinspin_prod(t2,g2,s2[ic2][ic1]);
	
	summ_the_trace_prod_spinspins(c,t1,t2);
      }
}

//Trace the product of gamma1 * spinspin1^dag * gamma2 * spinspin2,
void trace_g_sdag_g_s(complex *glb_c,dirac_matr *g1,colorspinspin *s1,dirac_matr *g2,colorspinspin *s2,const int ncontr)
{
  //Allocate a contiguous memory area where to store local node results
  complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
  for(int icontr=0;icontr<ncontr;icontr++)
    for(int glb_t=0;glb_t<glb_size[0];glb_t++) loc_c[icontr*glb_size[0]+glb_t][0]=loc_c[icontr*glb_size[0]+glb_t][1]=0;
  
  for(int icontr=0;icontr<ncontr;icontr++) 
    {
      if(debug_lvl>1) master_printf("Contraction %d/%d\n",icontr+1,ncontr);

      //Local loop
      for(int ivol=0;ivol<loc_vol;ivol++)
	{
	  int glb_t=glb_coord_of_loclx[ivol][0];
	  //Color loop
	  for(int ic=0;ic<3;ic++)
	    {
	      complex ctemp;
	      site_trace_g_sdag_g_s(ctemp,&(g1[icontr]),s1[ivol][ic],&(g2[icontr]),s2[ivol][ic]);
	      complex_summassign(loc_c[icontr*glb_size[0]+glb_t],ctemp);
	    }
	}
    }
  
  //Final reduction
  if(debug_lvl>1) master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(debug_lvl>1) master_printf("Reduction done\n");
  
  nissa_free(loc_c);
}
  
void sum_trace_g_sdag_g_s_times_trace_g_sdag_g_s(complex **glb_c, dirac_matr *g1L,colorspinspin *s1L, dirac_matr *g2L, colorspinspin *s2L, dirac_matr *g1R,colorspinspin *s1R, dirac_matr *g2R, colorspinspin *s2R,const int ncontr)
{
//Allocate a contguous memory area where to store local results
  complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
  for(int icontr=0;icontr<ncontr;icontr++)
    for(int glb_t=0;glb_t<glb_size[0];glb_t++) loc_c[icontr*glb_size[0]+glb_t][0]=loc_c[icontr*glb_size[0]+glb_t][1]=0;

  //Check if the global vector is contiguos
  complex *glb_c_buf=glb_c[0];
  int use_buf=0;
  for(int icontr=0;icontr<ncontr && use_buf==0;icontr++) use_buf=(glb_c[icontr]!=glb_c_buf+icontr*glb_size[0]);
  if(use_buf)
    {
      if(debug_lvl>1)
        {
          master_printf("Creating a temporary buffer for the contractions.\n");
          master_printf("Avoid it passing a 'glb_c' pointing to a contiguos memory area\n");
        }
      glb_c_buf=nissa_malloc("glb_c_buf",ncontr*glb_size[0],complex);
    }

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      if(debug_lvl>1) master_printf("Contraction %d/%d\n",icontr+1,ncontr);

      //Local loop
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
        {
          complex ctemp;
	  complex ctempL_color,ctempL;
	  complex ctempR_color,ctempR;
	  //Initialize to zero
	  ctempL[0]=0.;
	  ctempL[1]=0.;
	  ctempR[0]=0.;
	  ctempR[1]=0.;
          int glb_t=glb_coord_of_loclx[loc_site][0];
          //Color loops
          for(int icol=0;icol<3;icol++) {
		site_trace_g_sdag_g_s(ctempL_color,&(g1L[icontr]),s1L[loc_site][icol],&(g2L[icontr]),s2L[loc_site][icol]);
		complex_summ(ctempL,ctempL,ctempL_color);
         }
	 for(int icol=0;icol<3;icol++) {   
	      site_trace_g_sdag_g_s(ctempR_color,&(g1R[icontr]),s1R[loc_site][icol],&(g2R[icontr]),s2R[loc_site][icol]);
	      complex_summ(ctempR,ctempR,ctempR_color);
          }
	      safe_complex_prod(ctemp,ctempL,ctempR);
              complex_summ(loc_c[icontr*glb_size[0]+glb_t],loc_c[icontr*glb_size[0]+glb_t],ctemp);
        }
    }
    

  //Finale reduction
  if(debug_lvl>1) master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c_buf,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(debug_lvl>1) master_printf("Reduction done\n");

  //if a temporary buffer has been used, destory it after copyng data to the true one
  if(use_buf)
    {
      for(int icontr=0;icontr<ncontr;icontr++)
        for(int glb_t=0;glb_t<glb_size[0];glb_t++)
          {
            glb_c[icontr][glb_t][0]=glb_c_buf[icontr*glb_size[0]+glb_t][0];
            glb_c[icontr][glb_t][1]=glb_c_buf[icontr*glb_size[0]+glb_t][1];
          }
      nissa_free(glb_c_buf);
    }

  nissa_free(loc_c);
}


void trace_g_sdag_g_s_g_sdag_g_s(complex **glb_c, dirac_matr *g1L,colorspinspin *s1L, dirac_matr *g2L, colorspinspin *s2L, dirac_matr *g1R,colorspinspin *s1R, dirac_matr *g2R, colorspinspin *s2R,const int ncontr)
{
//Allocate a contguous memory area where to store local results
  complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
  for(int icontr=0;icontr<ncontr;icontr++)
    for(int glb_t=0;glb_t<glb_size[0];glb_t++) loc_c[icontr*glb_size[0]+glb_t][0]=loc_c[icontr*glb_size[0]+glb_t][1]=0;

  //Check if the global vector is contiguos
  complex *glb_c_buf=glb_c[0];
  int use_buf=0;
  for(int icontr=0;icontr<ncontr && use_buf==0;icontr++) use_buf=(glb_c[icontr]!=glb_c_buf+icontr*glb_size[0]);
  if(use_buf)
    {
      if(debug_lvl>1)
        {
          master_printf("Creating a temporary buffer for the contractions.\n");
          master_printf("Avoid it passing a 'glb_c' pointing to a contiguos memory area\n");
        }
      glb_c_buf=nissa_malloc("glb_c_buf",ncontr*glb_size[0],complex);
    }

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      if(debug_lvl>1) master_printf("Contraction %d/%d\n",icontr+1,ncontr);

      //Local loop
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
        {
          complex ctemp;
          int glb_t=glb_coord_of_loclx[loc_site][0];
          //Color loop
          for(int icol1=0;icol1<3;icol1++) for(int icol2=0;icol2<3;icol2++)
            {
//              site_trace_g_sdag_g_s_g_sdag_g_s(ctemp,&(g1L[icontr]),s1L[loc_site][icol1],&(g2L[icontr]),s2L[loc_site][icol2],&(g1R[icontr]),s1R[loc_site][icol2],&(g2R[icontr]),s2R[loc_site][icol1]);
		site_trace_g_sdag_g_s_g_sdag_g_s(ctemp,&(g1L[icontr]),s1L[loc_site][icol2],&(g2L[icontr]),s2R[loc_site][icol2],&(g1R[icontr]),s1R[loc_site][icol1],&(g2R[icontr]),s2L[loc_site][icol1]);

              complex_summ(loc_c[icontr*glb_size[0]+glb_t],loc_c[icontr*glb_size[0]+glb_t],ctemp);
            }
        }
    }

  //Final reduction
  if(debug_lvl>1) master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c_buf,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(debug_lvl>1) master_printf("Reduction done\n");

  //if a temporary buffer has been used, destory it after copyng data to the true one
  if(use_buf)
    {
      for(int icontr=0;icontr<ncontr;icontr++)
        for(int glb_t=0;glb_t<glb_size[0];glb_t++)
          {
            glb_c[icontr][glb_t][0]=glb_c_buf[icontr*glb_size[0]+glb_t][0];
            glb_c[icontr][glb_t][1]=glb_c_buf[icontr*glb_size[0]+glb_t][1];
          }
      nissa_free(glb_c_buf);
    }

  nissa_free(loc_c);
}

void trace_id_sdag_g_s_id_sdag_g_s(complex *glb_c,colorspinspin *s1L,dirac_matr *g2L,colorspinspin *s2L,colorspinspin *s1R,dirac_matr *g2R,colorspinspin *s2R,const int ncontr)
{
  //Allocate a contiguous memory area where to store local results
  complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
  memset(loc_c,0,sizeof(complex)*ncontr*glb_size[0]);

  //Local loop
  spinspin AR,AL;
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    {
      int glb_t=glb_coord_of_loclx[loc_site][0];
      for(int icol1=0;icol1<3;icol1++)
	for(int icol2=0;icol2<3;icol2++)
	  {
	    spinspin_spinspindag_prod(AL,s2L[loc_site][icol1],s1L[loc_site][icol2]);
	    spinspin_spinspindag_prod(AR,s2R[loc_site][icol2],s1R[loc_site][icol1]);
	    
	    for(int icontr=0;icontr<ncontr;icontr++)
	      {
                complex ctemp;
		spinspin ALg, ARg;
		
                spinspin_dirac_spinspin_prod(ALg,g2R+icontr,AL);
                spinspin_dirac_spinspin_prod(ARg,g2L+icontr,AR);
                trace_prod_spinspins(ctemp,ALg,ARg);
                complex_summassign(loc_c[icontr*glb_size[0]+glb_t],ctemp);
	      }
	  }
    }

  //Final reduction
  if(debug_lvl>1) master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(debug_lvl>1) master_printf("Reduction done\n");
  
  nissa_free(loc_c);
}

//Trace the product of gamma1 * su3spinspin1^dag * gamma2 * su3spinspin2,
void trace_g_ccss_dag_g_ccss(complex *glb_c,dirac_matr *g1,su3spinspin *s1,dirac_matr *g2,su3spinspin *s2,const int ncontr)
{
  //Allocate a contiguous memory area where to store local node results
  complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
  for(int icontr=0;icontr<ncontr;icontr++)
    for(int glb_t=0;glb_t<glb_size[0];glb_t++) loc_c[icontr*glb_size[0]+glb_t][0]=loc_c[icontr*glb_size[0]+glb_t][1]=0;
  
  for(int icontr=0;icontr<ncontr;icontr++) 
    {
      if(debug_lvl>1) master_printf("Contraction %d/%d\n",icontr+1,ncontr);

      //Local loop
      for(int ivol=0;ivol<loc_vol;ivol++)
        {
          int glb_t=glb_coord_of_loclx[ivol][0];
	  
	  complex ctemp;
	  site_trace_g_ccss_dag_g_ccss(ctemp,&(g1[icontr]),s1[ivol],&(g2[icontr]),s2[ivol]);
	  complex_summassign(loc_c[icontr*glb_size[0]+glb_t],ctemp);
        }
    }
  
  //Final reduction
  if(debug_lvl>1) master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(debug_lvl>1) master_printf("Reduction done\n");
  
  nissa_free(loc_c);
}

void sum_trace_id_sdag_g_s_times_trace_id_sdag_g_s(complex *glb_c,colorspinspin *s1L,dirac_matr *g2L,colorspinspin *s2L,colorspinspin *s1R,dirac_matr *g2R,colorspinspin *s2R,const int ncontr)
{
  //Allocate a contguous memory area where to store local results
  complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
  memset(loc_c,0,sizeof(complex)*ncontr*glb_size[0]);
  
  //Local loop
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    {
      int glb_t=glb_coord_of_loclx[loc_site][0];
      complex ctempL[ncontr];
      complex ctempR[ncontr];
      complex ctemp[ncontr];
      
      memset(ctempL,0,sizeof(complex)*ncontr);
      memset(ctempR,0,sizeof(complex)*ncontr);
      memset(ctemp,0,sizeof(complex)*ncontr);
      
      for(int icol=0;icol<3;icol++)
	{
	  spinspin AL,AR;
	  spinspin_spinspindag_prod(AL,s2L[loc_site][icol],s1L[loc_site][icol]);
	  spinspin_spinspindag_prod(AR,s2R[loc_site][icol],s1R[loc_site][icol]);
	  for(int icontr=0;icontr<ncontr;icontr++)
	    {
	      complex ctempL_color,ctempR_color;

	      trace_prod_dirac_spinspin(ctempL_color,g2R+icontr,AL);
	      trace_prod_dirac_spinspin(ctempR_color,&(g2L[icontr]),AR);

	      complex_summ(ctempL[icontr],ctempL[icontr],ctempL_color);
	      complex_summ(ctempR[icontr],ctempR[icontr],ctempR_color);
	    }
	}
      
      for(int icontr=0;icontr<ncontr;icontr++)
	complex_summ_the_prod(loc_c[icontr*glb_size[0]+glb_t],ctempL[icontr],ctempR[icontr]);
    }
  
  //Final reduction
  if(debug_lvl>1) master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(debug_lvl>1) master_printf("Reduction done\n");
  
  nissa_free(loc_c);
}

//print a single contraction to the file
void print_contraction_to_file(FILE *fout,int op1,int op2,complex *contr,int twall,const char *tag,double norm)
{
  if(rank==0)
    {
      fprintf(fout," # %s%s%s\n",tag,gtag[op2],gtag[op1]);
      for(int tempt=0;tempt<glb_size[0];tempt++)
	{
	  int t=tempt+twall;
	  if(t>=glb_size[0]) t-=glb_size[0];
	  
	  fprintf(fout,"%+016.16g\t%+016.16g\n",contr[t][0]*norm,contr[t][1]*norm);
	}
    }
}

//print all the passed contractions to the file
void print_contractions_to_file(FILE *fout,int ncontr,int *op1,int *op2,complex *contr,int twall,const char *tag,double norm)
{
  if(rank==0)
    for(int icontr=0;icontr<ncontr;icontr++)
      {
	fprintf(fout,"\n");
	print_contraction_to_file(fout,op1[icontr],op2[icontr],contr+icontr*glb_size[0],twall,tag,norm);
      }
}

//Rotate left and right by (1+-ig5)/sqrt(2)
//We distinguish for types of rotations, ir=rsink*2+rsource
//We divide the spinspin into 4 blocks, according to id<2
// -------------------------------
// | div |  0  |  1  |  2  |  3  |
// -------------------------------
// | 0 1 | + 1 | 1 + | 1 - | - 1 |
// | 2 3 | 1 - | - 1 | + 1 | 1 + |
// -------------------------------
// so just have to specify which is the block which rotate as +-i
void rotate_spinspin_to_physical_basis(spinspin s,int rsi,int rso)
{
  const int list_prb[4]={0,1,2,3},list_mrb[4]={3,2,1,0}; //plus and minus rotating blocks
  const int so_shft[4]={0,2,0,2},si_shft[4]={0,0,2,2};   //start of dirac indexes defining blocks
  
  int ir=rsi*2+rso,prb=list_prb[ir],mrb=list_mrb[ir];
  
  for(int dso=0;dso<2;dso++)
    for(int dsi=0;dsi<2;dsi++)
      {
	int pso=dso+so_shft[prb],psi=dsi+si_shft[prb];
	int mso=dso+so_shft[mrb],msi=dsi+si_shft[mrb];
	
	// rotation with +,-
	assign_complex_prod_i(s[pso][psi]);
	assign_complex_prod_minus_i(s[mso][msi]);
      }
}

void rotate_vol_colorspinspin_to_physical_basis(colorspinspin *s,int rsi,int rso)
{
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int ic=0;ic<3;ic++)
      rotate_spinspin_to_physical_basis(s[ivol][ic],rsi,rso);
}

void rotate_vol_su3spinspin_to_physical_basis(su3spinspin *s,int rsi,int rso)
{
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int ic1=0;ic1<3;ic1++)
      for(int ic2=0;ic2<3;ic2++)
	rotate_spinspin_to_physical_basis(s[ivol][ic1][ic2],rsi,rso);
}

//Calculate a lot of mesonic contractions
void lot_of_mesonic_contractions(complex *glb_contr,int **op,int ncontr,colorspinspin **S0,colorspinspin **S1,intpair npr,intpair *pr_combo,int npr_combo,int twall)
{
  //compute spatial volume
  int spat_loc_vol=loc_size[1]*loc_size[2]*loc_size[3];
  
  //compute the local buffer dimension
  int loc_buf_size=npr_combo*loc_size[0]*ncontr;
  
  //allocate local buffer
  complex *loc_contr=nissa_malloc("loc3pts",loc_buf_size,complex);
  memset(loc_contr,0,sizeof(complex)*loc_buf_size);
  
  //remap the operators listinto dirac matrixes, putting g5 to revert second prop
  //and compactifying the source/sink multiplication
  dirac_matr d[2][ncontr];
  int eq_to_mult[2][ncontr],nind_mult[2];
  //loop over source or sink and contraction
  for(int ss=0;ss<2;ss++)
    {
      //reset the independence and number of indep contractions
      for(int icontr=0;icontr<ncontr;icontr++) eq_to_mult[ss][icontr]=-1;
      nind_mult[ss]=0;
      
      //loop over all contractions
      for(int icontr=0;icontr<ncontr;icontr++)
	if(eq_to_mult[ss][icontr]==-1)
	  {
	    //if indep, save a new indep multiplications
	    if(ss==0) dirac_prod(&(d[ss][nind_mult[ss]]), &(base_gamma[op[ss][icontr]]),&(base_gamma[5]));
	    else      dirac_prod(&(d[ss][nind_mult[ss]]), &(base_gamma[5]),&(base_gamma[op[ss][icontr]]));
	    eq_to_mult[ss][icontr]=nind_mult[ss];
	    
	    //find equivalent ones
	    for(int jcontr=icontr+1;jcontr<ncontr;jcontr++)
	      if(eq_to_mult[ss][jcontr]==-1 && op[ss][jcontr]==op[ss][icontr]) eq_to_mult[ss][jcontr]=eq_to_mult[ss][icontr];
	    
	    //and increment the number of indep
	    nind_mult[ss]++;
	  }
    }
  
  //allocate room where to store the gamma times each prop at each point
  spinspin **ss[2];
  for(int iss=0;iss<2;iss++)
    {
      ss[iss]=nissa_malloc("cs_ptr",npr[iss],spinspin*);
      ss[iss][0]=nissa_malloc("ss",npr[iss]*nind_mult[iss],spinspin);
      for(int ip=1;ip<npr[iss];ip++) ss[iss][ip]=ss[iss][0]+nind_mult[iss]*ip;
    }
  
  //summ over all the spatial volume the contractions
  int ivol=0;
  //loop over t
  for(int t=0;t<loc_size[0];t++)
    for(int ispat=0;ispat<spat_loc_vol;ispat++)
      {
	for(int ic=0;ic<3;ic++)
	  {
	    //if(ic==0) master_printf("point %d over %d,time: %lg\n",ivol,loc_vol,take_time());
	    
	    //bufferize for each point all the propagators multiplied by all required gamma
	    for(int ip=0;ip<npr[0];ip++)
	      for(int im=0;im<nind_mult[0];im++)
		spinspin_dirac_spinspindag_prod(ss[0][ip][im],&(d[0][im]),S0[ip][ivol][ic]);
	    for(int ip=0;ip<npr[1];ip++)
	      for(int im=0;im<nind_mult[1];im++)
		spinspin_dirac_spinspin_prod_transp(ss[1][ip][im],&(d[1][im]),S1[ip][ivol][ic]);
	    
	    int offset=t;
	    //span all propagator combinations adding the full trace to the local contraction result
	    for(int ipr_combo=0;ipr_combo<npr_combo;ipr_combo++)
	      for(int icontr=0;icontr<ncontr;icontr++)
		{
		  complex *A=(complex*)ss[0][pr_combo[ipr_combo][0]][eq_to_mult[0][icontr]],*B=(complex*)ss[1][pr_combo[ipr_combo][1]][eq_to_mult[1][icontr]];
#ifdef BGP
		  bgp_complex cpu_out,cpu_A,cpu_B;
		  bgp_load_complex(cpu_out,loc_contr[offset]);
#pragma unroll(15)
		  for(int i=0;i<15;i++)
		    {
		      bgp_load_complex(cpu_A,A[i]);
		      bgp_load_complex(cpu_B,B[i]);
		      bgp_complex_summ_the_prod(cpu_out,cpu_out,__lfpd(A[i]),__lfpd(B[i]));
		      bgp_cache_touch_complex(A[i+1]);
		      bgp_cache_touch_complex(B[i+1]);
		    }
		  bgp_load_complex(cpu_A,A[15]);
		  bgp_load_complex(cpu_B,B[15]);
		  bgp_complex_summ_the_prod(cpu_out,cpu_out,__lfpd(A[15]),__lfpd(B[15]));
		  bgp_save_complex(loc_contr[offset],cpu_out);
#else
		  for(int i=0;i<16;i++)
		    complex_summ_the_prod(loc_contr[offset],A[i],B[i]);
#endif
		  
		  //increment the offset of the buffer
		  offset+=loc_size[0];
		}
	  }
	ivol++;
      }

  //find the number of propagator combinations to be hold on each x0=0 rank
  int nrank_x0=rank_tot/nrank_dir[0];
  int npr_combo_per_rank_x0=(int)ceil((double)npr_combo/nrank_x0);
  int nrank_x0_tbu=(int)ceil((double)npr_combo/npr_combo_per_rank_x0);
  //master_printf("Ncombo tot=%d, ncombo_per_rank_x0=%d, nrank to be used: %d\n",npr_combo,npr_combo_per_rank_x0,nrank_x0_tbu);
  int npr_combo_lst_rank_x0=npr_combo-npr_combo_per_rank_x0*(nrank_x0_tbu-1);
  
  //perform complanar ranks reduction
  //master_printf("Reducing...\n");
  for(int irank=0;irank<nrank_x0_tbu;irank++)
    {
      int npr_combo_cur_rank_x0=(irank<(nrank_x0_tbu-1)) ? npr_combo_per_rank_x0 : npr_combo_lst_rank_x0;
      //master_printf("Reducing on %d rank %d combo\n",irank,npr_combo_cur_rank_x0);
      MPI_Reduce(loc_contr+loc_size[0]*ncontr*npr_combo_per_rank_x0*irank,
		 glb_contr,loc_size[0]*ncontr*npr_combo_cur_rank_x0,
		 MPI_DOUBLE_COMPLEX,MPI_SUM,irank,plan_comm[0]);
    }
  //master_printf("Reduction done!\n");
  
  //finalize the operations on master rank of each line communicator of x0=0 plan
  if(plan_rank[0]<nrank_x0_tbu)
    {
      int npr_combo_cur_rank=(plan_rank[0]<nrank_x0_tbu-1) ? npr_combo_per_rank_x0 : npr_combo_lst_rank_x0;
      int npart_corr_cur_rank=npr_combo_cur_rank*ncontr;
      int loc_buf_size=npart_corr_cur_rank*loc_size[0];
      
      if(line_rank[0]!=0)
      	//on non-master node send data
	MPI_Send((void*)glb_contr,loc_buf_size,MPI_DOUBLE_COMPLEX,0,786+line_rank[0],line_comm[0]);
      else
	{
	  //on global master node open incoming communications
	  MPI_Request request[nrank_dir[0]];
	  for(int trank=1;trank<nrank_dir[0];trank++)
	    MPI_Irecv((void*)(glb_contr+loc_buf_size*trank),loc_buf_size,MPI_DOUBLE_COMPLEX,trank,786+trank,line_comm[0],&request[trank-1]);
	  
	  //prepare final reordering: transpose t with all the other indices and shift backward of twall
	  int *ord=nissa_malloc("reord",loc_buf_size*nrank_dir[0],int);
	  for(int trank=0;trank<nrank_dir[0];trank++)
	    for(int icorr=0;icorr<npart_corr_cur_rank;icorr++)
	      for(int loc_t=0;loc_t<loc_size[0];loc_t++)
		{
		  int sour_t=loc_t+loc_size[0]*trank;
		  int dest_t=(sour_t>=twall) ? sour_t-twall : sour_t+glb_size[0]-twall;
		  
		  int fr=loc_t+loc_size[0]*(icorr+npart_corr_cur_rank*trank);
		  int to=dest_t+glb_size[0]*icorr;
		  
		  //compute the swapped order index
		  ord[fr]=to;
		}
	  
	  //wait to finish receiving all data
	  //master_printf("Waiting to receive all data\n");
	  MPI_Waitall(nrank_dir[0]-1,request,MPI_STATUS_IGNORE);
	  //master_printf("All data received\n");
	  
	  //reorder the received buffer
	  //master_printf("Reordering data\n");
	  reorder_vector((char*)glb_contr,ord,loc_buf_size*nrank_dir[0],sizeof(complex));
	  nissa_free(ord);
	  //master_printf("Data reordered\n");
	}
    }
  
  //free all vectors
  nissa_free(ss[0][0]);
  nissa_free(ss[1][0]);
  
  nissa_free(ss[0]);
  nissa_free(ss[1]);

  nissa_free(loc_contr);
}

