#include "nissa.h"

/*
             SEQ	                  	 
            .....                        .....           
          ..     ..                    ..     ..         
         .         .                  .         .        
    op  X           X  g5   =    op  X           .  S1 
         .         .                  .         .        
          ..     ..                    ..     ..         
            .....                        .....           
              S0                                    
*/

//This function contract a source with a sequential spinor putting the passed list of operators
void contract_with_source(complex **corr,colorspinspin *S1,int *list_op,colorspinspin *source,int ncontr,int fprop,int rprop)
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
  trace_g_sdag_g_s(corr,t2,source,t1,S1,ncontr);
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  initNissa();

  //take initial time
  double tic;
  if(debug_lvl)
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
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);

  //Init the MPI grid 
  initGrid(T,L);

  //allocate the source and the prop
  colorspinspin *source=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  colorspinspin *S1=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  if((source==NULL || S1==NULL ) && rank==0)
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
  read_colorspinspin(source,path,NULL);

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

      if(rank==0 && debug_lvl) printf(" contr.%d %d\n",icontr,op[icontr]);
    }

  //Read the number of the props
  int nprop;
  read_str_int("NProp",&nprop);

  //Inizialization of output
  FILE *fout=NULL;
  if(rank==0) fout=fopen("two_points_check","w");
  int spat_vol=glb_size[1]*glb_size[1]*glb_size[1];

  //Loop over propagators
  double mass_spec,mass1;
  double theta_spec,theta1;
  int phys;
  int r1;
  for(int iprop=0;iprop<nprop;iprop++)
    {

      //Read the name, mass, theta and r for the list of the propagators
      read_str(path,1024);
      read_double(&mass_spec);
      read_double(&theta_spec);
      read_double(&mass1);
      read_double(&theta1);
      read_int(&phys);
      read_int(&r1);

      if(debug_lvl && rank==0)
	printf(" prop.%d %s, m_spec=%f th_spec=%f, m1=%f th1=%f phys=%d r1=%d\n",iprop,path,mass_spec,theta_spec,mass1,theta1,phys,r1);
      
      //Read the propagator
      read_colorspinspin(S1,path,NULL);
      
      if(rank==0) fprintf(fout," # m_spec=%f th_spec=%f m1=%f th1=%f r1=%d\n",mass_spec,theta_spec,mass1,theta1,r1);
      
      contract_with_source(contr,S1,op,source,ncontr,phys,r1);

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
  if(debug_lvl)
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
  
  free(S1);
  free(source);

  free(op);
  free(contr[0]);
  free(contr);

  closeNissa();

  return 0;
}
