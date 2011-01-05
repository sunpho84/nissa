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

//Trace the product of gamma1 * spinspin1^dag * gamma2 * spinspin2
void trace_g_sdag_g_s(complex **glb_c,dirac_matr *g1,colorspinspin *s1,dirac_matr *g2,colorspinspin *s2,const int ncontr)
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
