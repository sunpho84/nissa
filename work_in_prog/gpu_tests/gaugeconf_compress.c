#include "nissa.h"

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  char filename[1024];
  int T,L;
  read_str_int("L",&L);
  read_str_int("T",&T);
  read_str_str("Filename",filename,1024);
  
  close_input();

  //Init the MPI grid 
  init_grid(T,L);

  ///////////////////////////////////////////

  quad_su3 *conf=allocate_quad_su3(loc_vol,"conf");
  
  read_gauge_conf(conf,filename);  

  test_gaugeconf_compression(conf,loc_vol);
  
  close_nissa();

  return 0;
}
