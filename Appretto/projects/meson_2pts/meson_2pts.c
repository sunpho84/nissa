
//This program calculates a list of two-point functions contracting 
//all the propagators in the first list with all the propagators in
//the second list. The propagators of the second list are loaded one 
//by one. The first list is split into blocks, each of them as large 
//as possible.

//First list has to be the shortest, and should contain to be reverted
//The meson has quark content:
//            _
//            q(m1)q(m2)
//
//This is the reference scheme:

/*               +              
                S(m1)                           
                .....          
              ..     ..        
             .         .       
        op1  X           X  op2
             .         .       
              ..     ..        
                .....          
                S(m2)          
                                    
source |------>---->----->---->| sink

*/



#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

//Calculate the maximum number of allocable propagators
int compute_allocable_propagators(int nprop_list)
{
  //First of all check if there is enough room for two propagators.
  //This is the minimal requirement for the program to be able to work.
  colorspinspin *fuf;
  fuf=(colorspinspin*)malloc(2*sizeof(colorspinspin)*loc_vol);

  if(fuf==NULL)
    {
      if(rank==0)
	{
	  fprintf(stderr,"Error: not enough memory for two propagators\n");
	  fflush(stderr);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }
  else if(debug>1 && rank==0) printf("Ok there is enough memory to load two propagators\n");

  free(fuf);

  //Now determine the largest number of propagator of the first list (and one additional) loadable at once.
  //We are sure that we can allocate at least 2 props, so it will exit with at least nprop_max=1.
  int nprop_max=nprop_list + 1;
  do
    {
      nprop_max--;
      fuf=(colorspinspin*)malloc((nprop_max+1)*sizeof(colorspinspin)*loc_vol);
    }
  while(fuf==NULL);

  free(fuf);

  if(debug && rank==0)
    printf("Will allocate %d propagators from a list with %d propagators\n",nprop_max,nprop_list);

  return nprop_max;
}

//This function takes care to make the revert on the FIRST spinor, putting the needed gamma5
//It also applies the appropriate rotators to the physical basis if asked
void meson_two_points(complex **corr,int *list_op1,colorspinspin *s1,int *list_op2,colorspinspin *s2,int ncontr,int f1,int r1,int f2,int r2)
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

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_appretto();

  if(narg<2 && rank==0)
      {
	fprintf(stderr,"Use: %s input_file\n",arg[0]);
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD,1);
      }

  open_input(arg[1]);

  //Read the volume
  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));

  //Read the time location of the source
  int twall;
  read_str_int("TWall",&twall);

  //Read the number of propagators of the first list
  int nprop_list1;
  read_str_int("NPropFirstlist",&nprop_list1);
  if(rank==0) printf("Nprop of the first list: %d\n",nprop_list1);

  //Read the name, mass, theta and other flags for the first list
  char **base_filename1=(char**)malloc(sizeof(char*)*nprop_list1);
  for(int iprop1=0;iprop1<nprop_list1;iprop1++) base_filename1[iprop1]=(char*)malloc(1024);
  double  *mass_prop1=(double*)malloc(sizeof(double)*nprop_list1);
  double *theta_prop1=(double*)malloc(sizeof(double)*nprop_list1);
  int    * phys_prop1=   (int*)malloc(sizeof(int)   *nprop_list1);
  int        *r_prop1=   (int*)malloc(sizeof(int)   *nprop_list1);
  for(int iprop=0;iprop<nprop_list1;iprop++)
    {
      read_str(base_filename1[iprop],1024);
      read_double(&(mass_prop1[iprop]));
      read_double(&(theta_prop1[iprop]));
      read_int(&(phys_prop1[iprop]));
      read_int(&(r_prop1[iprop]));

      if(debug && rank==0)
	printf(" prop.%d %s, m=%f th=%f phys=%d r=%d\n",iprop,base_filename1[iprop],mass_prop1[iprop],theta_prop1[iprop],phys_prop1[iprop],r_prop1[iprop]);
    }
      
  //Read the number of propagators of the second list
  int nprop_list2;
  read_str_int("NPropSecondlist",&nprop_list2);
  if(rank==0) printf("Nprop of the second list: %d\n",nprop_list2);

  //Read the name, mass, theta and other flags for the second list
  char **base_filename2=(char**)malloc(sizeof(char*)*nprop_list2);
  for(int iprop2=0;iprop2<nprop_list2;iprop2++) base_filename2[iprop2]=(char*)malloc(1024);
  double  *mass_prop2=(double*)malloc(sizeof(double)*nprop_list2);
  double *theta_prop2=(double*)malloc(sizeof(double)*nprop_list2);
  int    * phys_prop2=   (int*)malloc(sizeof(int)   *nprop_list2);
  int        *r_prop2=   (int*)malloc(sizeof(int)   *nprop_list2);
  for(int iprop=0;iprop<nprop_list2;iprop++)
    {
      read_str(base_filename2[iprop],1024);
      read_double(&(mass_prop2[iprop]));
      read_double(&(theta_prop2[iprop]));
      read_int(&(phys_prop2[iprop]));
      read_int(&(r_prop2[iprop]));

      if(debug && rank==0)
	printf(" prop.%d %s, m=%f th=%f phys=%d r=%d\n",iprop,base_filename2[iprop],mass_prop2[iprop],theta_prop2[iprop],phys_prop2[iprop],r_prop2[iprop]);
    }
      
  //Read the number of contractions
  int ncontr;
  read_str_int("NContr",&ncontr);
  if(rank==0) printf("Number of contractions: %d\n",ncontr);

  //Initialize the list of correlations and the list of operators
  //contiguous allocation
  complex **contr=(complex**)malloc(sizeof(complex*)*ncontr);
  contr[0]=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]); 
  int *op1=(int*)malloc(sizeof(int)*ncontr);
  int *op2=(int*)malloc(sizeof(int)*ncontr);
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      contr[icontr]=contr[0]+icontr*glb_size[0];

      //Read the operator pairs
      read_int(&(op1[icontr]));
      read_int(&(op2[icontr]));

      if(rank==0 && debug) printf(" contr.%d %d %d\n",icontr,op1[icontr],op2[icontr]);
    }
  
  //Read the output filename
  char outfile[1024];
  read_str_str("Output",outfile,1024);

  close_input();

  //Init the MPI grid 
  init_grid();
  
  //Calculate the number of blocks for the first list
  int nprop_per_block=compute_allocable_propagators(nprop_list1);
  int nblocks=nprop_list1/nprop_per_block;
  if(nprop_list1>nblocks*nprop_per_block) nblocks++;

  //allocate the spinors
  colorspinspin **spinor1=(colorspinspin**)malloc(sizeof(colorspinspin*)*nprop_per_block);
  for(int iprop1=0;iprop1<nprop_per_block;iprop1++) spinor1[iprop1]=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  colorspinspin *spinor2=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);

  ///////////////////////////////////////////
  
  //take initial time
  double tic;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }

  FILE *fout=NULL;
  if(rank==0)
    {
      fout=fopen(outfile,"w");
      if(fout==NULL)
	{
	  fprintf(stderr,"Couldn't open the file: %s",outfile);
	  fflush(stderr);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }
  int spat_vol=glb_size[1]*glb_size[1]*glb_size[1];

  //Loop over the blocks of the first list
  for(int iblock=0;iblock<nblocks;iblock++)
    {
      int iblock_first=iblock*nprop_per_block;
      int iblock_last=min_int((iblock+1)*nprop_per_block,nprop_list1);
      int iblock_length=iblock_last-iblock_first;

      if(rank==0 && debug) printf("Block %d/%d length: %d\n",iblock,nblocks,iblock_length);

      //now read the whole first block
      for(int iprop1=0;iprop1<iblock_length;iprop1++)
      {
	int counter=iblock_first+iprop1;
	
	if(debug>1 && rank==0) printf("Going to read propagator %d/%d: %s\n",iprop1,iblock_length,base_filename1[counter]);
	read_colorspinspin(base_filename1[counter],spinor1[iprop1]);
      }

      //now loop over the second popagator
      for(int iprop2=0;iprop2<nprop_list2;iprop2++)
	{
	  //read the second propagator one by one
	  colorspinspin *spinor2_ptr; //This will point to spinor2 if the prop. is not in the first list
	  //check if the file is already loaded
	  spinor2_ptr=spinor2;
	  for(int iprop1=0;iprop1<iblock_length;iprop1++)
	    {
	      int counter=iblock_first+iprop1;
	      if(strcmp(base_filename1[counter],base_filename2[iprop2])==0)
		{
		  spinor2_ptr=spinor1[iprop1];
		  if(debug && rank==0) printf("Propagator %s found in the position %d of the first list\n",base_filename2[iprop2],counter);
		}
	    }
	  //if not found in the first list, load it
	  if(spinor2_ptr==spinor2)
	    {
	      if(debug>1 && rank==0) printf("Going to read propagator %d: %s\n",iprop2,base_filename2[iprop2]);
	      read_colorspinspin(base_filename2[iprop2],spinor2);
	    }
	  
	  //Calculate all the two points between spinor 1 and spinor2
	  for(int iprop1=0;iprop1<iblock_length;iprop1++)
	    {
	      int counter=iblock_first+iprop1;

	      if(rank==0)
		fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d\n",mass_prop1[counter],theta_prop1[counter],r_prop1[counter],mass_prop2[iprop2],theta_prop2[iprop2],r_prop2[iprop2]);

	      meson_two_points(contr,op1,spinor1[iprop1],op2,spinor2_ptr,ncontr,phys_prop1[counter],r_prop1[counter],phys_prop2[iprop2],r_prop2[iprop2]);

	      if(rank==0)
		{
		  fprintf(fout,"\n");
		  
		  for(int icontr=0;icontr<ncontr;icontr++)
		    {
		      fprintf(fout," # %s%s\n",gtag[op2[icontr]],gtag[op1[icontr]]);
		      fprintf(fout,"\n");
		      for(int tempt=0;tempt<glb_size[0];tempt++)
			{
			  int t=tempt+twall;
			  if(t>=glb_size[0]) t-=glb_size[0];
			  
			  fprintf(fout,"%+016.16g\t%+016.16f\n",contr[icontr][t][0]/spat_vol,contr[icontr][t][1]/spat_vol);
			}
		      fprintf(fout,"\n");
		    }
		}
	      if(rank==0)
		{
		  fprintf(fout,"\n");
		  if(debug>1) fflush(fout);
		}
	    }
	}
    }

  //take final time
  double tac;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tac=MPI_Wtime();
      if(rank==0) printf("Time elapsed in doing %d contractions: %f s\n",nprop_list1*nprop_list2*ncontr,tac-tic);
    }

  ///////////////////////////////////////////

  if(rank==0)
    {
      printf("\nEverything ok, exiting!\n");
      fclose(fout);
    }

  free(mass_prop2);
  free(theta_prop2);
  free(phys_prop2);
  free(r_prop2);
  free(spinor2);
  for(int iprop2=0;iprop2<nprop_list2;iprop2++) free(base_filename2[iprop2]);
  free(base_filename2);

  free(mass_prop1);
  free(theta_prop1);
  free(phys_prop1);
  free(r_prop1);
  for(int iprop1=0;iprop1<nprop_per_block;iprop1++) free(spinor1[iprop1]);
  free(spinor1);
  for(int iprop1=0;iprop1<nprop_list1;iprop1++) free(base_filename1[iprop1]);
  free(base_filename1);

  free(contr[0]);
  free(contr);

  free(op1);
  free(op2);

  close_appretto();

  return 0;
}
