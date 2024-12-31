#include "nissa.hpp"

using namespace nissa;

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<4) CRASH("use: %s L T file_out",arg[0]);
  
  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  
  //Init the MPI grid
  init_grid(T,L);
  
  //crete and write
  quad_su3 *conf=nissa_malloc("conf",locVol,quad_su3);
  generate_cold_lx_conf(conf);
  NISSA_LOC_VOL_LOOP(site)
    for(int dir=0;dir<4;dir++)
      su3_prodassign_double(conf[site][dir],glbCoordOfLoclx[site][dir]);
  
  write_ildg_gauge_conf(arg[3],conf,64);
  
  nissa_free(conf);
  
  ///////////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
