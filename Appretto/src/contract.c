#pragma once

//Local trace of the product of gamma1 * spinspin1^dag * gamma2 * spinspin2
void site_trace_g_s_g_sdag(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2)
{
  spinspin t1,t2;
  
  spinspin_dirac_spinspin_prod(t1,g1,s1);
  spinspin_dirac_spinspindag_prod(t2,g2,s2);

  trace_prod_spinspins(c,t1,t2);
}

//Trace the product of gamma1 * spinspin1^dag * gamma2 * spinspin2
void trace_g_s_g_sdag(complex *glb_c,dirac_matr *g1,colorspinspin *s1,dirac_matr *g2,colorspinspin *s2,const int ncontr)
{
  complex loc_c[ncontr][glb_size[0]];

  for(int icontr=0;icontr<ncontr;icontr++) 
    {
      //Reset the output
      for(int t=0;t<glb_size[0];t++)
	loc_c[icontr][t][0]=loc_c[icontr][t][1]=0;
      
      //Local loop
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
	{
	  complex ctemp;
	  int glb_t=glb_coord[loc_site][0];
	  //Color loop
	  for(int icol=0;icol<3;icol++)
	    {
	      site_trace_g_s_g_sdag(ctemp,&(g1[icontr]),s1[loc_site][icol],&(g2[icontr]),s2[loc_site][icol]);
	      complex_summ(loc_c[icontr][glb_t],loc_c[icontr][glb_t],ctemp);
	    }
	}
    }
  
  //Finale reduction
  MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
}


