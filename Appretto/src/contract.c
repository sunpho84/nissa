#pragma once

//Local trace of the product of gamma1 * spinspin1^dag * gamma2 * spinspin2
void site_trace_g_sdag_g_s(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2)
{
  spinspin t1,t2;
  
  spinspin_dirac_spinspindag_prod(t1,g1,s1);
  spinspin_dirac_spinspin_prod(t2,g2,s2);

  trace_prod_spinspins(c,t1,t2);
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
      if(debug>1 && rank==0) printf("Contrction %d/%d\n",icontr+1,ncontr);
      
      //Local loop
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
	{
	  complex ctemp;
	  int glb_t=glb_coord[loc_site][0];
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


