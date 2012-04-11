#include "nissa.h"

int main(int narg,char **arg)
{
  double precision;
  char conf_in_path[1024];
  char conf_out_path[1024];
  
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //basic mpi initialization
  init_nissa();
  
  //Init the MPI grid 
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
  open_input(arg[1]);
  
  read_str_str("InGaugePath",conf_in_path,1024);
  read_str_double("Precision",&precision);
  read_str_str("OutGaugePath",conf_out_path,1024);
  
  close_input();
  
  ///////////////////////////////////////////
  
  quad_su3 *conf=nissa_malloc("Conf",loc_vol+bord_vol,quad_su3);
  quad_su3 *fix_conf=nissa_malloc("Conf2",loc_vol+bord_vol,quad_su3);
  
  read_ildg_gauge_conf(conf,conf_in_path);
  communicate_lx_quad_su3_borders(conf);
  
  landau_gauge_fix(fix_conf,conf,precision);
  
  write_gauge_conf(conf_out_path,fix_conf);
  
  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(conf));
  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(fix_conf));

  ///////////////////////////////////////////

  nissa_free(conf);
  nissa_free(fix_conf);
  
  close_nissa();

  return 0;
}
