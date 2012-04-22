#include "nissa.h"

int main(int narg,char **arg)
{
  char base_filename1[1024],base_filename2[1024];

  //basic mpi initialization
  init_nissa();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  //read lattice size
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //init the grid
  init_grid(T,L);

  //read the propagators filename
  read_str_str("BaseFilename1",base_filename1,1024);
  read_str_str("BaseFilename2",base_filename2,1024);
  
  //read the number of contractions
  int ncontr;
  read_str_int("NContr",&ncontr);
  
  //read the list of operators pair
  //available operators [0-15] = S, V1,V2,V3,V4, P, A1,A2,A3,A4, T1,T2,T3, B1,B2,B3 
  expect_str("ListSourceSinkOpPairs");
  dirac_matr op1[ncontr],op2[ncontr];
  int iop1[ncontr],iop2[ncontr];
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      read_int(&(iop1[icontr]));
      read_int(&(iop2[icontr]));
      op1[icontr]=base_gamma[iop1[icontr]];
      op2[icontr]=base_gamma[iop2[icontr]];
    }

  close_input();

  ///////////////////////////////////////////
  
  //allocate vectors for contractions
  complex **contr=(complex**)malloc(sizeof(complex*)*ncontr);
  for(int icontr=0;icontr<ncontr;icontr++)
    contr[icontr]=(complex*)malloc(sizeof(complex)*glb_size[0]);
  
  colorspinspin *spinor1=allocate_colorspinspin(loc_vol,"spinor1"); //perform automatically also allocation test
  colorspinspin *spinor2=allocate_colorspinspin(loc_vol,"spinor2"); //the string is used as ref. in case of error

  double load_time=-take_time();
  char *filename_ending=NULL; //no suffix
  read_colorspinspin(spinor1,base_filename1,filename_ending);
  read_colorspinspin(spinor2,base_filename2,filename_ending);
  load_time+=take_time(); 

  //perform the contraction
  double contr_time=-take_time();
  trace_g_sdag_g_s(contr,op1,spinor1,op2,spinor2,ncontr);
  contr_time+=take_time();

  if(rank==0)
    {
      printf("Time elapsed in loading: %g s\n",load_time);
      printf("Time elapsed in contracting: %g s\n",contr_time);
    }

  //this is needed only if the source is a time wall
  int spat_vol=glb_size[1]*glb_size[2]*glb_size[3];
  if(rank==0)
    for(int icontr=0;icontr<ncontr;icontr++)
      {
	printf("Contraction: %s-%s\n",gtag[iop1[icontr]],gtag[iop2[icontr]]);
	for(int t=0;t<glb_size[0];t++)
	  printf("t=%d %+16.16g %+16.16g\n",t,contr[icontr][t][0]/spat_vol,contr[icontr][t][1]/spat_vol);
	printf("\n");
      }

  check_free(spinor1);
  check_free(spinor2);

  for(int icontr=0;icontr<ncontr;icontr++)
    check_free(contr[icontr]);
  
  check_free(contr);
  
  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
