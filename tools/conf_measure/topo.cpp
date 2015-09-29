#include "nissa.hpp"

using namespace nissa;

int L,T;

THREADABLE_FUNCTION_1ARG(unitarize_conf_max, quad_su3*,conf)
{
  GET_THREAD_ID();
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int idir=0;idir<4;idir++)
      {
	su3 t;
	su3_unitarize_orthonormalizing(t,conf[ivol][idir]);
	su3_copy(conf[ivol][idir],t);
      }
  set_borders_invalid(conf);
}
THREADABLE_FUNCTION_END


void in_main(int narg,char **arg)
{
  if(narg<2) crash("use: %s input",arg[0]);
  
  //open input file
  open_input(arg[1]);
  
  //init the grid
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
  //read in and out conf path
  char conf_path[1024];
  read_str_str("ConfPath",conf_path,1024);
  
  //read topo pars
  top_meas_pars_t top_meas_pars;
  read_top_meas_pars(top_meas_pars,true);
  
  //////////////////////////// read the conf /////////////////////////////
  
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  //read the conf and write plaquette
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  read_ildg_gauge_conf(conf,conf_path,&mess);
  unitarize_conf_max(conf);
  
  measure_topology_lx_conf(top_meas_pars,conf,0,0);
  
  nissa_free(conf);
  ILDG_message_free_all(&mess);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
