#pragma once

void trace_prod_dirac_spinspin(complex c,dirac_matr *a,spinspin b)
{
  c[0]=c[1]=0;
  for(int id=0;id<4;id++)
    complex_summ_the_prod(c,a->entr[id],b[id][a->pos[id]]);
}

void spinspin_spinspindag_prod(spinspin out,spinspin a,spinspin b)
{
  //This is the line on the matrix
  for(int id1=0;id1<4;id1++)
    for(int id2=0;id2<4;id2++)
      {
        out[id1][id2][0]=out[id1][id2][1]=0;
        for(int id=0;id<4;id++) complex_summ_the_conj2_prod(out[id1][id2],a[id1][id],b[id2][id]);
      }
}

void trace_id_sdag_g_s_id_sdag_g_s(complex *glb_c,colorspinspin *s1L,dirac_matr *g2L,colorspinspin *s2L,colorspinspin *s1R,dirac_matr *g2R,colorspinspin *s2R,const int ncontr)
{
  //Allocate a contiguous memory area where to store local results
  complex *loc_c=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]);
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
		
                spinspin_dirac_spinspin_prod(ALg,&(g2R[icontr]),AL);
                spinspin_dirac_spinspin_prod(ARg,&(g2L[icontr]),AR);
                trace_prod_spinspins(ctemp,ALg,ARg);
                complex_summassign(loc_c[icontr*glb_size[0]+glb_t],ctemp);
	      }
	  }
    }

  //Final reduction
  if(debug>1 && rank==0) printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(debug>1 && rank==0) printf("Reduction done\n");
  
  check_free(loc_c);
}

void sum_trace_id_sdag_g_s_times_trace_id_sdag_g_s(complex *glb_c,colorspinspin *s1L,dirac_matr *g2L,colorspinspin *s2L,colorspinspin *s1R,dirac_matr *g2R,colorspinspin *s2R,const int ncontr)
{
  //Allocate a contguous memory area where to store local results
  complex *loc_c=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]);
  memset(loc_c,0,sizeof(complex)*ncontr*glb_size[0]);

  //Local loop
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    {
      int glb_t=glb_coord_of_loclx[loc_site][0];
      complex ctempL[ncontr];
      complex ctempR[ncontr];
      complex ctemp[ncontr];
      
      for(int icontr=0;icontr<ncontr;icontr++)
	for(int reim=0;reim<2;reim++)
	  ctempL[icontr][reim]=ctempR[icontr][reim]=ctemp[icontr][reim]=0;
      
      for(int icol=0;icol<3;icol++)
	{
	  spinspin AL;
	  spinspin_spinspindag_prod(AL,s2L[loc_site][icol],s1L[loc_site][icol]);
	  for(int icontr=0;icontr<ncontr;icontr++)
	    {
	      complex ctempL_color;
	      trace_prod_dirac_spinspin(ctempL_color,&(g2R[icontr]),AL);
	      complex_summ(ctempL[icontr],ctempL[icontr],ctempL_color);
	    }
	}
      for(int icol=0;icol<3;icol++)
	{
	  spinspin AR;
	  spinspin_spinspindag_prod(AR,s2R[loc_site][icol],s1R[loc_site][icol]);
	  
	  for(int icontr=0;icontr<ncontr;icontr++)
	    {
	      complex ctempR_color;
	      trace_prod_dirac_spinspin(ctempR_color,&(g2L[icontr]),AR);
	      complex_summ(ctempR[icontr],ctempR[icontr],ctempR_color);
	    }
	}
      for(int icontr=0;icontr<ncontr;icontr++)
	{
	  safe_complex_prod(ctemp[icontr],ctempL[icontr],ctempR[icontr]);
	  complex_summ(loc_c[icontr*glb_size[0]+glb_t],loc_c[icontr*glb_size[0]+glb_t],ctemp[icontr]);
      }
    }
  
  //Final reduction
  if(debug>1 && rank==0) printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(debug>1 && rank==0) printf("Reduction done\n");
  
  check_free(loc_c);
}
