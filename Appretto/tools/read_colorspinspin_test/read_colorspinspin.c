#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

int main(int narg,char **arg)
{
  char base_filename[1024];

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
  read_str_str("BaseFilename",base_filename,1024);

  close_input();

  //Init the MPI grid 
  init_grid();
  
  ///////////////////////////////////////////

  colorspinspin *spinore=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);

  read_colorspinspin(spinore,base_filename);

  if(rank==0)
    for(int loc_site=0;loc_site<loc_vol;loc_site++)
      for(int ic=0;ic<3;ic++)
	{
	  for(int id_sink=0;id_sink<4;id_sink++)
	    {
	      for(int id_source=0;id_source<4;id_source++)
		printf("%g,%g\t",spinore[loc_site][ic][id_source][id_sink][0],spinore[loc_site][ic][id_source][id_sink][0]);
	      printf("\n");
	    }
	      printf("\n");
	}
	      
  free(spinore);

  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
