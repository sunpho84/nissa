#pragma once

#include <stdio.h>

#include "global.c"

//Initilialise the random number generator
void init_random(int seed)
{
  const int im1=2147483563,ia1=40014;
  const int iq1=53668,ir1=12211;
  int j,k;

  //Generate the true seed
  if(rank==0)
    {
      srand(seed);
      seed=rand();
      if(debug>1) printf("New seed generated, %d\n",seed);
    }
  MPI_Bcast(&seed,1,MPI_INT,0,MPI_COMM_WORLD);

  if(rank==0 && loc_vol==0)
    {
      fprintf(stderr,"Error, grid not initalized\n");
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  
  ran2_idum=(int*)malloc(loc_vol*sizeof(int));
  ran2_idum2=(int*)malloc(loc_vol*sizeof(int));
  ran2_iy=(int*)malloc(loc_vol*sizeof(int));
  ran2_iv=(int**)malloc(loc_vol*sizeof(int*));
  for(int loc_ind=0;loc_ind<loc_vol;loc_ind++)
    {
      ran2_iv[loc_ind]=(int*)malloc(ran2_ntab*sizeof(int));
      
      //initialization
      ran2_idum[loc_ind]=seed+glblx_of_loclx[loc_ind];

      ran2_idum[loc_ind]=max_int(ran2_idum[loc_ind]+1,1);
      ran2_idum2[loc_ind]=ran2_idum[loc_ind];
      for(j=ran2_ntab+7;j>=0;j--)
        {
          k=ran2_idum[loc_ind]/iq1;
          ran2_idum[loc_ind]=ia1*(ran2_idum[loc_ind]-k*iq1)-k*ir1;
          if(ran2_idum[loc_ind]<0) ran2_idum[loc_ind]+=im1;
          if(j<ran2_ntab) ran2_iv[loc_ind][j]=ran2_idum[loc_ind];
        }
      ran2_iy[loc_ind]=ran2_iv[loc_ind][0];
    }

  if(rank==0 && debug) printf("Random generator initialized\n");
}

//Close random
void close_random()
{
  free(ran2_idum);
  free(ran2_idum2);
  free(ran2_iy);
  for(int loc_ind=0;loc_ind<loc_vol;loc_ind++) free(ran2_iv[loc_ind]);
  free(ran2_iv);
}

//standard ran2 from numerical recipes
double ran2(int loc_ind)
{
  const int im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014,ia2=40692;
  const int iq1=53668,iq2=52774,ir1=12211,ir2=3791,ndiv=1+imm1/ran2_ntab;
  const double am=1.0/im1,eps=1.2e-7,rnmx=1-eps;
  int j,k;
  double out;

  k=ran2_idum[loc_ind]/iq1;
  ran2_idum[loc_ind]=ia1*(ran2_idum[loc_ind]-k*iq1)-k*ir1;
  if(ran2_idum[loc_ind]<0) ran2_idum[loc_ind]+=im1;
  
  k=ran2_idum2[loc_ind]/iq2;
  ran2_idum2[loc_ind]=ia2*(ran2_idum2[loc_ind]-k*iq2)-k*ir2;
  if(ran2_idum2[loc_ind]<0) ran2_idum2[loc_ind]+=im2;
      
  j=ran2_iy[loc_ind]/ndiv;
  ran2_iy[loc_ind]=ran2_iv[loc_ind][j]-ran2_idum2[loc_ind];
  ran2_iv[loc_ind][j]=ran2_idum[loc_ind];
  if(ran2_iy[loc_ind]<0) ran2_iy[loc_ind]+=imm1;

  out=min_double(am*ran2_iy[loc_ind],rnmx);

  return out;
}

//return +-1
int pm_one(int loc_site)
{
  double r=ran2(loc_site);
  if(r>0.5) return 1;
  else return -1;
}

