#include "nissa.h"

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();

  //initialize the program
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //read the input
  open_input(arg[1]);

  //Read the volume
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L); 
  //Kappa
  double kappa;
  read_str_double("Kappa",&kappa);
  //Mass
  double mass;
  read_str_double("Mass",&mass);
  //Residue
  double residue;
  read_str_double("Residue",&residue);
  //Path of conf
  char gauge_conf_path[1024];
  read_str_str("GaugeConfPath",gauge_conf_path,1024);
  
  close_input();
  
  //read conf
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  read_ildg_gauge_conf(conf,gauge_conf_path);
  
  //prepare the total source
  spincolor *source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  memset(source,0,loc_vol*sizeof(spincolor));
  if(rank==0) source[0][0][0][0]=1;
  set_borders_invalid(source);
  
  //allocate solution
  spincolor *sol=nissa_malloc("sol",loc_vol+bord_vol,spincolor);
  
  inv_tmQ2_cg_128(sol,NULL,conf,kappa,mass,10000,residue,source);
  
  nissa_free(sol);
  nissa_free(source);
  nissa_free(conf);
  
  close_nissa();
  
  return 0;
}
