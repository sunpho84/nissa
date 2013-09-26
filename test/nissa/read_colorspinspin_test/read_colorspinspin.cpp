#include "nissa.h"

int main(int narg,char **arg)
{
  char base_filename[1024];

  //basic mpi initialization
  init_nissa();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  //read the sizes
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L);
  
  read_str_str("BaseFilename",base_filename,1024);

  close_input();

  ///////////////////////////////////////////

  colorspinspin *spinore=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);

  read_colorspinspin(spinore,base_filename,NULL);

  if(rank==0)
    NISSA_LOC_VOL_LOOP(ivol)
      for(int ic=0;ic<3;ic++)
	{
	  printf(" # %d %d\n",ivol,ic);
	  for(int id_sink=0;id_sink<4;id_sink++)
	    {
	      for(int id_source=0;id_source<4;id_source++)
		printf("%g,%g\t",spinore[ivol][ic][id_sink][id_source][0],spinore[ivol][ic][id_sink][id_source][1]);
	      printf("\n");
	    }
	      printf("\n");
	}
	      
  free(spinore);

  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
