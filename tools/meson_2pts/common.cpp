#pragma once

#include "nissa.h"

//Calculate the maximum number of allocable propagators
//First of all check if there is enough room for the configuration,
//then for the two/three propagators.
//This is the minimal requirement for the program to be able to work.
int compute_allocable_propagators(int nprop_list,int nch_contr,int nmin_req)
{
  quad_su3 *temp_conf=NULL;
  as2t_su3 *temp_clov=NULL;
  if(nch_contr>0)
    {
      nmin_req++;
      
      temp_conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol+bord_vol+edge_vol));
      if(temp_conf==NULL && rank>0)
	{
	  fprintf(stderr,"Unable to allocate the space for the gauge configuration!\n");
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

      temp_clov=(as2t_su3*)malloc(sizeof(as2t_su3)*loc_vol);
      if(temp_clov==NULL && rank>0)
	{
	  fprintf(stderr,"Unable to allocate the space for the P_munu term!\n");
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }

  colorspinspin *fuf=NULL;

  fuf=(colorspinspin*)malloc(nmin_req*sizeof(colorspinspin)*loc_vol);

  if(fuf==NULL)
    {
      if(rank==0)
	{
	  fprintf(stderr,"Error: not enough memory for %d propagators\n",nmin_req);
	  fflush(stderr);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }
  else if(debug_lvl>1 && rank==0) printf("Ok there is enough memory to load %d propagators\n",nmin_req);

  free(fuf);

  //Now determine the largest number of propagator of the first list (and one additional) loadable at once.
  //We are sure that we can allocate at least nmin_req props, so it will exit with at least nprop_max=1.
  int nprop_max=nprop_list+nmin_req-1;
  do
    {
      nprop_max--;
      fuf=(colorspinspin*)malloc((nprop_max+1)*sizeof(colorspinspin)*loc_vol);
    }
  while(fuf==NULL);

  free(fuf);

  if(debug_lvl && rank==0)
    printf("Will allocate %d propagators from a list with %d propagators\n",nprop_max,nprop_list);

  if(nch_contr>0)
    {
      free(temp_conf);
      free(temp_clov);
    }

  return nprop_max;
}

//This function takes care to make the revert on the FIRST spinor, putting the needed gamma5
//It also applies the appropriate rotators to the physical basis if asked
void meson_two_points(complex *corr,int *list_op1,colorspinspin *s1,int *list_op2,colorspinspin *s2,int ncontr,int f1,int r1,int f2,int r2)
{
  //Temporary vectors for the internal gamma
  dirac_matr t1[ncontr],t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      
      //Put the two gamma5 needed for the revert of the first spinor
      dirac_prod(&(t1[icontr]), &(base_gamma[list_op1[icontr]]),&(base_gamma[5]));
      dirac_prod(&(t2[icontr]), &(base_gamma[5]),&(base_gamma[list_op2[icontr]]));
      
      //Remind that D- rotates as 1+ig5, but D-^-1 rotates as 1-ig5,
      //moreover (D^-1)^dagger rotates again as 1+ig5 (pweee!!!)

      //f1 < 1: do rotation for a quark propagator
      //f1 = 1: do nothing (physical basis already)
      //f1 > 1: do rotation for a (charged) sequential propagator

      if(f1<1)
	switch(r1)
	  {
	  case 0: //This is (D-^-1)^dag
	    dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pplus);
	    dirac_prod(&(t2[icontr]), &Pplus,&(t2[icontr]));
	    break;
	  case 1: //This is (D+^-1)^dag
	    dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pminus);
	    dirac_prod(&(t2[icontr]), &Pminus,&(t2[icontr]));
	    break;
	  }

      if(f1>1)
        switch(r1)
          {
          case 0: //This is (D-^-1)^dag
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pplus);
            dirac_prod(&(t2[icontr]), &Pminus,&(t2[icontr]));
            break;
          case 1: //This is (D+^-1)^dag
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pminus);
            dirac_prod(&(t2[icontr]), &Pplus,&(t2[icontr]));
            break;
          }

      if(f2<1)
	switch(r2)
	  {
	  case 0: //This is D-^-1
	    dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pminus);
	    dirac_prod(&(t1[icontr]), &Pminus,&(t1[icontr]));
	    break;
	  case 1: //This is D+^-1
	    dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pplus);
	    dirac_prod(&(t1[icontr]), &Pplus,&(t1[icontr]));
	    break;
	  }

     if(f2>1)
        switch(r2)
          {
          case 0: //This is D-^-1
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pminus);
            dirac_prod(&(t1[icontr]), &Pplus,&(t1[icontr]));
            break;
          case 1: //This is D+^-1
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pplus);
            dirac_prod(&(t1[icontr]), &Pminus,&(t1[icontr]));
            break;
          }

    }
  
  //Call for the routine which does the real contraction
  trace_g_sdag_g_s(corr,t1,s1,t2,s2,ncontr);
}

