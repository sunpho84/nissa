#include "nissa.hpp"

using namespace nissa;


void in_main(int narg,char **arg)
{
  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  init_grid(T,L);
  
  start_loc_rnd_gen(atoi(arg[4]));
  
  ///////////////////////////////////////////
  
  quad_su3 *conf=nissa_malloc("Conf",locVol+bord_vol,quad_su3);
  quad_su3 *fix_conf=nissa_malloc("Conf2",locVol+bord_vol,quad_su3);
  
  read_ildg_gauge_conf(conf,arg[3]);
  communicate_lx_quad_su3_borders(conf);
  
  perform_random_gauge_transform(fix_conf,conf);
  
  write_ildg_gauge_conf(arg[5],fix_conf,64);
  
  master_printf("plaq: %16.16lg\n",global_plaquette_lx_conf(conf));
  master_printf("plaq: %16.16lg\n",global_plaquette_lx_conf(fix_conf));
  
  ///////////////////////////////////////////
  
  nissa_free(conf);
  nissa_free(fix_conf);
  
}

int main(int narg,char **arg)
{
  if(narg<6) crash("Use: %s L T path_in seed path_out",arg[0]);
  
  //basic mpi initialization
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
