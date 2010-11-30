#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

/*
            	 
           .....           
         ..     ..         
        .         .        
   op  X           .  S
        .         .        
         ..     ..         
           .....           
                      
*/

//This function contract a source with a propagator putting the passed list of operators
void contract_with_source(complex **corr,colorspinspin *prop,int *list_op,colorspinspin *source,int ncontr,int fprop,int rprop)
{
  //Temporary vector for the internal matrices
  dirac_matr t1[ncontr],t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      
      //Init the second 
      t1[icontr]=base_gamma[list_op[icontr]];
      t2[icontr]=base_gamma[0];
      
      //Remind that D- rotates as 1+ig5, but D-^-1 rotates as 1-ig5,

      //f1 < 1: do rotation for a quark propagator
      //f1 = 1: do nothing (physical basis already)
      //f1 > 1: do rotation for a (charged) sequential propagator

      if(fprop<1)
	switch(rprop)
	  {
	  case 0: //This is D-^-1
	    dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pminus);
	    dirac_prod(&(t2[icontr]), &Pminus,&(t2[icontr]));
	    break;
	  case 1: //This is D+^-1
	    dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pplus);
	    dirac_prod(&(t2[icontr]), &Pplus,&(t2[icontr]));
	    break;
	  }

      if(fprop>1)
        switch(rprop)
          {
          case 0: //This is D-^-1 on the sequential. t1 is on the sink
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pminus);
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pplus);
            break;
          case 1: //This is D+^-1 on the sequential. t1 is on the sink
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pplus);
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pminus);
            break;
          }

    }

  //Call for the routine which does the real contraction
  trace_g_sdag_g_s(corr,t2,source,t1,prop,ncontr);
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_appretto();

  //take initial time
  double tic;
  if(debug)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      tic=MPI_Wtime();
    }

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

  //Init the MPI grid 
  init_grid();

  //allocate the source and the prop
  colorspinspin *source=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  colorspinspin *prop=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  if((source==NULL || prop==NULL ) && rank==0)
    {
      fprintf(stderr,"Error! Not enoug memory to allocate source and spinor!\n");
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  //Read the time location of the source
  int twall;
  read_str_int("TWall",&twall);

  //Read the source
  char path[1024];
  read_str_str("Source",path,1024);
  read_colorspinspin(source,path);

  //Read the number of contractions
  int ncontr;
  read_str_int("NContr",&ncontr);
  if(rank==0) printf("Number of contractions: %d\n",ncontr);

  //Initialize the list of correlations and the list of operators
  complex **contr=(complex**)malloc(sizeof(complex*)*ncontr);
  contr[0]=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]);
  int *op=(int*)malloc(sizeof(int)*ncontr);
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      contr[icontr]=contr[0]+icontr*glb_size[0];
      //Read the operator pairs
      read_int(&(op[icontr]));

      if(rank==0 && debug) printf(" contr.%d %d\n",icontr,op[icontr]);
    }

  //Inizialization of output
  char outfile[1024];
  FILE *fout=NULL;
  read_str_str("Output",outfile,1024);
  if(rank==0) fout=fopen(outfile,"w");
  int spat_vol=glb_size[1]*glb_size[1]*glb_size[1];

  //Read the number of the props
  int nprop;
  read_str_int("NProp",&nprop);

  //Loop over propagators
  int phys,r;
  char tag[1024];
  for(int iprop=0;iprop<nprop;iprop++)
    {
      //Read the name, mass, theta and r for the list of the propagators
      read_str(path,1024);
      read_str(tag,1024);
      read_int(&phys);
      read_int(&r);

      if(debug && rank==0)
	printf(" prop.%d %s tagged as '%s', phys=%d r=%d\n",iprop,path,tag,phys,r);
      
      //Read the propagator
      read_colorspinspin(prop,path);
      
      if(rank==0) fprintf(fout," # %s r=%d\n",tag,r);
      
      contract_with_source(contr,prop,op,source,ncontr,phys,r);

      if(rank==0)
	{
	  fprintf(fout,"\n");
	  
	  for(int icontr=0;icontr<ncontr;icontr++)
	    {
	      fprintf(fout," # %s\n\n",gtag[op[icontr]]);
	      for(int tempt=0;tempt<glb_size[0];tempt++)
		{
		  int t=tempt+twall;
		  if(t>=glb_size[0]) t-=glb_size[0];
		  
		  fprintf(fout,"%+016.16g\t%+016.16g\n",contr[icontr][t][0]/spat_vol,contr[icontr][t][1]/spat_vol);
		}
	      fprintf(fout,"\n");
	    }
	  fprintf(fout,"\n");
	}
    }

  close_input();

  //take final time
  double tac;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tac=MPI_Wtime();
      if(rank==0) printf("Time elapsed in doing %d contractions: %g s\n",nprop*ncontr,tac-tic);
    }

  ///////////////////////////////////////////
  
  if(rank==0)
    {
      printf("\nEverything ok, exiting!\n");
      fclose(fout);
    }
  
  free(prop);
  free(source);

  free(op);
  free(contr[0]);
  free(contr);

  close_appretto();

  return 0;
}
