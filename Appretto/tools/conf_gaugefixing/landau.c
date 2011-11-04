#include "appretto.h"

int main(int narg,char **arg)
{
  char conf_in_path[1024];
  double precision;

  //basic mpi initialization
  init_appretto();

  if(narg<2) crash("Use: %s input_file",arg[0]);

  open_input(arg[1]);

  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  read_str_str("GaugeConf",conf_in_path,1024);
  read_str_double("Precision",&precision);
  
  close_input();

  //Init the MPI grid 
  init_grid();

  ///////////////////////////////////////////

  quad_su3 *conf=appretto_malloc("Conf",loc_vol+loc_bord,quad_su3);
  quad_su3 *fix_conf=appretto_malloc("Conf2",loc_vol+loc_bord,quad_su3);
  
  read_gauge_conf(conf,conf_in_path);
  communicate_gauge_borders(conf);
  
  landau_gauge_fix(fix_conf,conf,precision);
  
  master_printf("plaq: %.18g\n",global_plaquette(conf));
  master_printf("plaq: %.18g\n",global_plaquette(fix_conf));

  ///////////////////////////////////////////

  appretto_free(conf);
  appretto_free(fix_conf);
  
  close_appretto();

  return 0;
}
