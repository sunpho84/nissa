#include "nissa.hpp"

using namespace nissa;

void in_main(int narg,char **arg)
{
  open_input(arg[1]);
  
  //Init the MPI grid
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
  char conf_in_path[1024];
  read_str_str("InGaugePath",conf_in_path,1024);
  double precision;
  read_str_double("Precision",&precision);
  char conf_out_path[1024];
  read_str_str("OutGaugePath",conf_out_path,1024);
  
  close_input();
  
  start_loc_rnd_gen(1000);
  
  ///////////////////////////////////////////
  
  quad_su3 *conf=nissa_malloc("Conf",locVol+bord_vol,quad_su3);
  quad_su3 *fix_conf=nissa_malloc("Conf2",locVol+bord_vol,quad_su3);
  
  read_ildg_gauge_conf(conf,conf_in_path);
  communicate_lx_quad_su3_borders(conf);
  
  //set pars
  LC_gauge_fixing_pars_t pars;
  pars.gauge=LC_gauge_fixing_pars_t::Landau;
  pars.target_precision=precision;
  
  crash("reimplement");
  //Landau_or_Coulomb_gauge_fix(fix_conf,&pars,conf);
  
  write_ildg_gauge_conf(conf_out_path,fix_conf,64);
  
  master_printf("plaq before: %16.16lg\n",global_plaquette_lx_conf(conf));
  master_printf("plaq after: %16.16lg\n",global_plaquette_lx_conf(fix_conf));
  
  ///////////////////////////////////////////
  
  nissa_free(conf);
  nissa_free(fix_conf);
  
}

int main(int narg,char **arg)
{
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //basic mpi initialization
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
