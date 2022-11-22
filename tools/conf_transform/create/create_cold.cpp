#include "nissa.hpp"

using namespace nissa;

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<4) crash("use: %s L T file_out",arg[0]);
  
  const int L=atoi(arg[1]);
  const int T=atoi(arg[2]);
  
  //Init the MPI grid
  init_grid(T,L);
  
  //crete and write
  LxField<quad_su3> conf("conf");
  generate_cold_lx_conf(conf);
  write_ildg_gauge_conf(arg[3],conf,64);
  
  ///////////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
