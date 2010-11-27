#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

int main(int narg,char **arg)
{
  char base_filename1[1024],base_filename2[1024];

  //basic mpi initialization
  init_appretto();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  read_str_str("BaseFilename1",base_filename1,1024);
  read_str_str("BaseFilename2",base_filename2,1024);

  close_input();

  //Init the MPI grid 
  init_grid();
  
  ///////////////////////////////////////////

  colorspinspin *spinore1=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  colorspinspin *spinore2=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);

  complex *contr=(complex*)malloc(sizeof(complex)*glb_size[0]);

  read_colorspinspin(base_filename1,spinore1);
  read_colorspinspin(base_filename2,spinore2);
  
  //take initial time
  double tic;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }
  trace_g_sdag_g_s(contr,&(base_gamma[0]),spinore1,&(base_gamma[0]),spinore2,1);
  //take final time
  double tac;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tac=MPI_Wtime();
      if(rank==0) printf("Time elapsed in contracting: %g d\n",tac-tic);
    }
  
  int spat_vol=glb_size[1]*glb_size[1]*glb_size[1];
  if(rank==0) for(int t=0;t<glb_size[0];t++) printf("%d %g %g",t,contr[t][0]/spat_vol,contr[t][1]/spat_vol);

  free(spinore1);
  free(spinore2);

  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
