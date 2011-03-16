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

  quad_su3 *origi_conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"or_conf");
  quad_su3 *smear_conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"sm_conf");
  
  read_local_gauge_conf(origi_conf,filename);
      
  ape_smearing(smear_conf,origi_conf,10,0.5);

  su3_print(origi_conf[0][1]);
  su3_print(smear_conf[0][1]);

  double origi_plaq=global_plaquette(origi_conf);
  double smear_plaq=global_plaquette(smear_conf);
  
  printf("%g %g\n",origi_plaq,smear_plaq);
  
  ///////////////////////////////////////////
  
  close_appretto();

  return 0;
}
