#pragma once

//Local trace of the product of gamma1 * spinspin1^dag * gamma2 * spinspin2
void site_trace_g_sdag_g_s(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2)
{
  spinspin t1,t2;
  
  spinspin_dirac_spinspindag_prod(t1,g1,s1);
  spinspin_dirac_spinspin_prod(t2,g2,s2);

  trace_prod_spinspins(c,t1,t2);
}
void site_trace_g_sdag_g_s_g_sdag_g_s (complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2, dirac_matr *g3,spinspin s3,dirac_matr *g4,spinspin s4)
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

//Trace the product of gamma1 * spinspin1^dag * gamma2 * spinspin2,
//this version is quite inefficient and kept here for pedagogical reasons
void trace_g_sdag_g_s_old(complex **glb_c,dirac_matr *g1,colorspinspin *s1,dirac_matr *g2,colorspinspin *s2,const int ncontr)
{
  //Allocate a contguous memory area where to store local results
  complex *loc_c=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]);
  for(int icontr=0;icontr<ncontr;icontr++)
    for(int glb_t=0;glb_t<glb_size[0];glb_t++) loc_c[icontr*glb_size[0]+glb_t][0]=loc_c[icontr*glb_size[0]+glb_t][1]=0;

  //Check if the global vector is contiguos
  complex *glb_c_buf=glb_c[0];
  int use_buf=0;
  for(int icontr=0;icontr<ncontr && use_buf==0;icontr++) use_buf=(glb_c[icontr]!=glb_c_buf+icontr*glb_size[0]);
  if(use_buf)
    {
      if(debug>1 && rank==0)
	{
	  printf("Creating a temporary buffer for the contractions.\n");
	  printf("Avoid it passing a 'glb_c' pointing to a contiguos memory area\n");
	}
      glb_c_buf=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]);
    }

  for(int icontr=0;icontr<ncontr;icontr++) 
    {
      if(debug>1 && rank==0) printf("Contraction %d/%d\n",icontr+1,ncontr);
      
      //Local loop
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
	{
	  complex ctemp;
	  int glb_t=glb_coord_of_loclx[loc_site][0];
	  //Color loop
	  for(int icol=0;icol<3;icol++)
	    {
	      site_trace_g_sdag_g_s(ctemp,&(g1[icontr]),s1[loc_site][icol],&(g2[icontr]),s2[loc_site][icol]);
	      complex_summ(loc_c[icontr*glb_size[0]+glb_t],loc_c[icontr*glb_size[0]+glb_t],ctemp);
	    }
	}
    }

  //Final reduction
  if(debug>1 && rank==0) printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c_buf,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(debug>1 && rank==0) printf("Reduction done\n");
  
  //if a temporary buffer has been used, destory it after copyng data to the true one
  if(use_buf)
    {
      for(int icontr=0;icontr<ncontr;icontr++)
	for(int glb_t=0;glb_t<glb_size[0];glb_t++)
	  {
	    glb_c[icontr][glb_t][0]=glb_c_buf[icontr*glb_size[0]+glb_t][0];
	    glb_c[icontr][glb_t][1]=glb_c_buf[icontr*glb_size[0]+glb_t][1];
	  }
      free(glb_c_buf);
    }

  free(loc_c);
}

//Trace the product of gamma1 * spinspin1^dag * gamma2 * spinspin2
//Performed this way:
// \sum_{vol,a,b,c,d} g1_ab S_bc^+ g2_cd S_da = \sum_{a,c} g1_a g2_c \sum_vol S_cb(a)^* S_d(c)a = 
// \sum_i X_i Y_i where i=4*a+c, X_i = g1_a g2_c and Y_i = \sum_vol S_cb(a)^* S_d(c)a
void trace_g_sdag_g_s(complex **glb_c,dirac_matr *g1,colorspinspin *s1,dirac_matr *g2,colorspinspin *s2,const int ncontr)
{
  //Allocate X_i
  complex X[loc_size[0]][16];
  //Contiguous memory area where to store local results
  complex loc_c[ncontr][glb_size[0]];
  memset(loc_c,0,glb_size[0]*ncontr*sizeof(complex));
  int T0=loc_size[0]*proc_coord[0];

  for(int icontr=0;icontr<ncontr;icontr++) 
    {
      if(debug>1 && rank==0) printf("Contraction %d/%d\n",icontr+1,ncontr);
      
      memset(X,0,loc_size[0]*16*sizeof(complex));
      
      //Local loop
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
	{
	  int loc_t=loc_coord_of_loclx[loc_site][0];

	  int i=0;
	  for(int a=0;a<4;a++)
	    for(int c=0;c<4;c++)
	      {
#ifdef BGP
		static _Complex double A,B,C;
		bgp_load_complex(A,X[loc_t][i]);
#endif
		//Color loop
		for(int icol=0;icol<3;icol++)
		  {
#ifdef BGP
		    bgp_load_complex(B,s1[loc_site][icol][c][g1[icontr].pos[a]]);
		    bgp_load_complex(C,s2[loc_site][icol][g2[icontr].pos[c]][a]);
		    bgp_complex_summ_the_conj1_prod(A,A,B,C);
#else
		    complex_summ_the_conj1_prod(X[loc_t][i],
						s1[loc_site][icol][c][g1[icontr].pos[a]],
						s2[loc_site][icol][g2[icontr].pos[c]][a]);
#endif
		  }
#ifdef BGP
		bgp_save_complex(X[loc_t][i],A);
#endif
		i++;
	      }
	    }
      
      int i=0;
      for(int a=0;a<4;a++)
	for(int c=0;c<4;c++)
	  {	
#ifdef BGP
	    static _Complex double Y,L,G1,G2,XI;
	    bgp_load_complex(G1,g1[icontr].entr[a]);
	    bgp_load_complex(G2,g2[icontr].entr[c]);
	    bgp_complex_prod(Y,G1,G2);
#else
	    complex Y;
	    unsafe_complex_prod(Y,g1[icontr].entr[a],g2[icontr].entr[c]);
#endif
	    for(int loc_t=0;loc_t<loc_size[0];loc_t++)
	      {
#ifdef BGP
		bgp_load_complex(L,loc_c[icontr][loc_t+T0]);
		bgp_load_complex(XI,X[loc_t][i]);
		bgp_complex_summ_the_prod(L,L,XI,Y);
		bgp_save_complex(loc_c[icontr][loc_t+T0],L);
#else
		complex_summ_the_prod(loc_c[icontr][loc_t+T0],X[loc_t][i],Y);
#endif
	      }
#ifdef BGP

#endif
	    i++;
	  }
    }

  //Check if the global vector is contiguous
  complex *glb_c_buf=glb_c[0];
  int use_buf=0;
  for(int icontr=0;icontr<ncontr && use_buf==0;icontr++) use_buf=(glb_c[icontr]!=glb_c_buf+icontr*glb_size[0]);
  if(use_buf)
    {
      if(debug>1 && rank==0)
	{
	  printf("Creating a temporary buffer for the contractions.\n");
	  printf("Avoid it passing a 'glb_c' pointing to a contiguos memory area\n");
	}
      glb_c_buf=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]);
    }

  //Final reduction
  if(debug>1 && rank==0) printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c_buf,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(debug>1 && rank==0) printf("Reduction done\n");
  
  //if a temporary buffer has been used, destory it after copying data to the true one
  if(use_buf)
    {
      for(int icontr=0;icontr<ncontr;icontr++)
	for(int glb_t=0;glb_t<glb_size[0];glb_t++)
	  {
	    glb_c[icontr][glb_t][0]=glb_c_buf[icontr*glb_size[0]+glb_t][0];
	    glb_c[icontr][glb_t][1]=glb_c_buf[icontr*glb_size[0]+glb_t][1];
	  }
      free(glb_c_buf);
    }
}

void sum_trace_g_sdag_g_s_times_trace_g_sdag_g_s(complex **glb_c, dirac_matr *g1L,colorspinspin *s1L, dirac_matr *g2L, colorspinspin *s2L, dirac_matr *g1R,colorspinspin *s1R, dirac_matr *g2R, colorspinspin *s2R,const int ncontr)
{
//Allocate a contguous memory area where to store local results
  complex *loc_c=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]);
  for(int icontr=0;icontr<ncontr;icontr++)
    for(int glb_t=0;glb_t<glb_size[0];glb_t++) loc_c[icontr*glb_size[0]+glb_t][0]=loc_c[icontr*glb_size[0]+glb_t][1]=0;

  //Check if the global vector is contiguos
  complex *glb_c_buf=glb_c[0];
  int use_buf=0;
  for(int icontr=0;icontr<ncontr && use_buf==0;icontr++) use_buf=(glb_c[icontr]!=glb_c_buf+icontr*glb_size[0]);
  if(use_buf)
    {
      if(debug>1 && rank==0)
        {
          printf("Creating a temporary buffer for the contractions.\n");
          printf("Avoid it passing a 'glb_c' pointing to a contiguos memory area\n");
        }
      glb_c_buf=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]);
    }

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      if(debug>1 && rank==0) printf("Contraction %d/%d\n",icontr+1,ncontr);

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
  if(debug>1 && rank==0) printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c_buf,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(debug>1 && rank==0) printf("Reduction done\n");

  //if a temporary buffer has been used, destory it after copyng data to the true one
  if(use_buf)
    {
      for(int icontr=0;icontr<ncontr;icontr++)
        for(int glb_t=0;glb_t<glb_size[0];glb_t++)
          {
            glb_c[icontr][glb_t][0]=glb_c_buf[icontr*glb_size[0]+glb_t][0];
            glb_c[icontr][glb_t][1]=glb_c_buf[icontr*glb_size[0]+glb_t][1];
          }
      free(glb_c_buf);
    }

  free(loc_c);
}


void trace_g_sdag_g_s_g_sdag_g_s(complex **glb_c, dirac_matr *g1L,colorspinspin *s1L, dirac_matr *g2L, colorspinspin *s2L, dirac_matr *g1R,colorspinspin *s1R, dirac_matr *g2R, colorspinspin *s2R,const int ncontr)
{
//Allocate a contguous memory area where to store local results
  complex *loc_c=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]);
  for(int icontr=0;icontr<ncontr;icontr++)
    for(int glb_t=0;glb_t<glb_size[0];glb_t++) loc_c[icontr*glb_size[0]+glb_t][0]=loc_c[icontr*glb_size[0]+glb_t][1]=0;

  //Check if the global vector is contiguos
  complex *glb_c_buf=glb_c[0];
  int use_buf=0;
  for(int icontr=0;icontr<ncontr && use_buf==0;icontr++) use_buf=(glb_c[icontr]!=glb_c_buf+icontr*glb_size[0]);
  if(use_buf)
    {
      if(debug>1 && rank==0)
        {
          printf("Creating a temporary buffer for the contractions.\n");
          printf("Avoid it passing a 'glb_c' pointing to a contiguos memory area\n");
        }
      glb_c_buf=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]);
    }

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      if(debug>1 && rank==0) printf("Contraction %d/%d\n",icontr+1,ncontr);

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

  //Finale reduction
  if(debug>1 && rank==0) printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c_buf,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(debug>1 && rank==0) printf("Reduction done\n");

  //if a temporary buffer has been used, destory it after copyng data to the true one
  if(use_buf)
    {
      for(int icontr=0;icontr<ncontr;icontr++)
        for(int glb_t=0;glb_t<glb_size[0];glb_t++)
          {
            glb_c[icontr][glb_t][0]=glb_c_buf[icontr*glb_size[0]+glb_t][0];
            glb_c[icontr][glb_t][1]=glb_c_buf[icontr*glb_size[0]+glb_t][1];
          }
      free(glb_c_buf);
    }

  free(loc_c);
}

//print all the passed contractions to the file
void print_contractions_to_file(FILE *fout,int ncontr,int *op1,int *op2,complex **contr,int twall,const char *tag)
{
  int spat_vol=glb_size[1]*glb_size[2]*glb_size[3];
  
  if(rank==0)
    for(int icontr=0;icontr<ncontr;icontr++)
      {
	fprintf(fout,"\n");
	fprintf(fout," # %s%s%s\n\n",tag,gtag[op2[icontr]],gtag[op1[icontr]]);
	for(int tempt=0;tempt<glb_size[0];tempt++)
	  {
	    int t=tempt+twall;
	    if(t>=glb_size[0]) t-=glb_size[0];
	    
	    fprintf(fout,"%+016.16g\t%+016.16g\n",contr[icontr][t][0]/spat_vol,contr[icontr][t][1]/spat_vol);
	  }
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
  const int so_shft[4]={0,2,0,2},si_shft[4]={0,0,2,2}; //shift of dirac indexes for blocks
  
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
