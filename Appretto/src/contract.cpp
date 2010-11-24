#pragma once

//Local trace of the product of gamma1 * spinspin1^dag * gamma2 * spinspin2
void site_trace_g_sdag_g_s(complex c,dirac_matr &g1,spinspin &s1,dirac_matr &g2,spinspin &s2)
{
  spinspin t1,t2;
  
  spinspin_dirac_spinspindag_prod(t1,g1,s1);
  spinspin_dirac_spinspin_prod(t2,g2,s2);

  trace_prod_spinspins(c,t1,t2);
}

//Trace of the product of gamma1 * spinspin1^dag * gamma2 * spinspin2
void trace_g_sdag_g_s(complex *glb_c,dirac_matr &g1,colorspinspin *s1,dirac_matr &g2,colorspinspin *s2)
{
  complex loc_c[glb_size[0]];
  for(int t=0;t<glb_size[0];t++) loc_c[t][0]=loc_c[t][1]=0;

  //Local loop
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    {
      complex ctemp;
      int glb_t=glb_coord[loc_site][0];
      //Color loop
      for(int icol=0;icol<3;icol++)
	{
	  site_trace_g_sdag_g_s(ctemp,g1,s1[loc_site][icol],g2,s2[loc_site][icol]);
	  complex_summ(loc_c[glb_t],loc_c[glb_t],ctemp);
	}
    }
  
  //Finale reduction
  MPI_Reduce(loc_c,glb_c,2*glb_size[0],MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
}

//Frontend to the gamma calculation
void contract_two_spinspins(complex *corr,dirac_matr &g1,colorspinspin *s1,dirac_matr &g2,colorspinspin *s2,int f1,int f2)
{
  dirac_matr t1,t2;

  //Put the two gamma5 needed for the revert of the first spinor
  dirac_prod(t1, g1,base_gamma[5]);
  dirac_prod(t2, base_gamma[5],g2);
  
  //Put the rotators to the physical basis these can be avoided by
  //putting a number different from 0 or 1 in f1 and f2

  //Remember that D- rotate as 1+ig5, but D-^-1 rotate ad 1-ig5,
  //morover (D^-1)^dagger rotate again as 1+ig5 (pweee!!!)

  if(f1==0) //This is (D-^-1)^dagger
    {
      dirac_prod(t1, t1,Pplus);
      dirac_prod(t2, Pplus,t2);
    }
  else      //This is (D+^-1)^dagger
    {
      dirac_prod(t1, t1,Pminus);
      dirac_prod(t2, Pminus,t2);
    }

  if(f2==0)  //This is D-^-1
    {
      dirac_prod(t2, t2,Pminus);
      dirac_prod(t1, Pminus,t1);
    }
  else       //This is D+^-1
    {
      dirac_prod(t2, t2,Pplus);
      dirac_prod(t1, Pplus,t1);
    }

  trace_g_sdag_g_s(corr,t1,s1,t2,s2);
}
