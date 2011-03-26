#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

int main(int narg,char **arg)
{
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
  
  close_input();

  //Init the MPI grid 
  init_grid();

  ///////////////////////////////////////////

  quad_su3 *origi_conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"or_conf");
  quad_su3 *smear_conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"sm_conf");
  
  read_local_gauge_conf(origi_conf,"/home/francesco/Prace/nissa/Appretto/test/gaugeconf_load/conf.0048");
      
  ape_smearing(smear_conf,origi_conf,0.4,7);
  
  su3_print(origi_conf[0][1]);
  su3_print(smear_conf[0][1]);

  double origi_plaq=global_plaquette(origi_conf);
  double smear_plaq=global_plaquette(smear_conf);
  
  printf("%g %g\n",origi_plaq,smear_plaq);
  
  ///////////////////////////////////////////
  
  spincolor *s=allocate_spincolor(loc_vol,"s");
  spincolor *t=allocate_spincolor(loc_vol,"t");
  read_spincolor(s,"/home/francesco/Prace/Programs/src/ahmidas-rw/test/point_src.48");
  jacobi_smearing(t,s,smear_conf,0.5,5);
  int l=loclx_of_coord_list(0,1,0,2);
  
  for(int d=0;d<4;d++)
    for(int c=0;c<3;c++)
      printf("%g %g\n",t[l][d][c][0],t[l][d][c][1]);
  
  close_appretto();

  return 0;
}
