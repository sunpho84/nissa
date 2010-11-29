#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

int main(int narg,char **arg)
{
  char filename[1024];

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
  read_str_str("Filename",filename,1024);

  close_input();

  //Init the MPI grid 
  init_grid();
  
  ///////////////////////////////////////////

  quad_su3 *conf=(quad_su3*)malloc(sizeof(quad_su3)*loc_vol);

  read_gauge_conf(conf,filename);

  free(conf);
  
  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
