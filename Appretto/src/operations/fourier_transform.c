#pragma once

//produce the table of the momentum
void Momentum(int **iP,double *bc,double *P2,double *SinP2,double **P,double **SinP,double *SinP4,int nmom)
{
  int lP[4],imom=0;

  for(lP[0]=iP[0][0];lP[0]<=iP[0][1];lP[0]++)
    for(lP[1]=iP[1][0];lP[1]<=iP[1][1];lP[1]++)
      for(lP[2]=iP[2][0];lP[2]<=iP[2][1];lP[2]++)
	for(lP[3]=iP[3][0];lP[3]<=iP[3][1];lP[3]++)
	  {
	    P2[imom]=SinP2[imom]=SinP4[imom]=0;
	    for(int idir=0;idir<4;idir++)
	      {
		P[idir][imom]=2*M_PI*(lP[idir]+bc[idir]*0.5)/glb_size[idir];
		P2[imom]+=pow(P[idir][imom],2);
		SinP[idir][imom]=sin(P[idir][imom]);
		SinP2[imom]+=pow(SinP[idir][imom],2);
		SinP4[imom]+=pow(SinP[idir][imom],4);
	      }
	    imom++;
	  }

  if(rank==0 && imom!=nmom)
    {
      fprintf(stderr,"imom != nmom\n");
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
}

//perform the fourier transform momentum per momentum
void spincolor_FT(spincolor *S,spincolor *FT,double *theta,int **iP,int nmom)
{
  double *P2=appretto_malloc("P2",nmom,double);
  double *SinP2=appretto_malloc("SinP2",nmom,double);
  double *SinP4=appretto_malloc("SinP4",nmom,double);

  double *P[4];
  double *SinP[4];

  for(int idir=0;idir<4;idir++)
    {
      P[idir]=appretto_malloc("P[idir]",nmom,double);
      SinP[idir]=appretto_malloc("SinP[idir]",nmom,double);
    }
  
  Momentum(iP,theta,P2,SinP2,P,SinP,SinP4,nmom); 

  //local FT
  spincolor *FT_loc=appretto_malloc("FT_loc",nmom,spincolor);
  for(int imom=0;imom<nmom;imom++) spincolor_put_to_zero(FT_loc[imom]); 
	
  double lP[4];
  for(int imom=0;imom<nmom;imom++)
    {	
      for(int idir=0;idir<4;idir++) lP[idir]=P[idir][imom];

      for(int ivol=0;ivol<loc_vol;ivol++)
	{
	  double arg=0;
	  for(int idir=0;idir<4;idir++) arg+=lP[idir]*glb_coord_of_loclx[ivol][idir];
	  complex eipx={cos(arg),-sin(arg)};

	  spincolor_summ_the_prod_complex(FT_loc[imom],S[ivol],eipx);			
	}
    }

  MPI_Reduce(FT_loc,FT,24*nmom,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  
  for(int idir=0;idir<4;idir++)
    {
      appretto_free(P[idir]);
      appretto_free(SinP[idir]);
    }
  appretto_free(P2);appretto_free(SinP2);appretto_free(SinP4);
}
