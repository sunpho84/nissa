#include "nissa.hpp"

using namespace nissa;

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<5) crash("use: %s L T seed file_out",arg[0]);

  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  int seed=atoi(arg[3]);

  //Init the MPI grid 
  init_grid(T,L);
  start_loc_rnd_gen(seed);

  //crete and write
  quad_su3 *conf=nissa_malloc("conf",locVol,quad_su3);
  generate_hot_lx_conf(conf);
  write_ildg_gauge_conf(arg[4],conf,64);
  
  nissa_free(conf);

  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
