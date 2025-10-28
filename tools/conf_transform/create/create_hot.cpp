#include "nissa.hpp"

using namespace nissa;

int main(int narg,char **arg)
{
  //basic mpi initialization
  initNissa(narg,arg);
  
  if(narg<5) CRASH("use: %s L T seed file_out",arg[0]);
  
  const int L=atoi(arg[1]);
  const int T=atoi(arg[2]);
  const int seed=atoi(arg[3]);
  
  //Init the MPI grid
  initGrid(T,L);
  start_loc_rnd_gen(seed);
  
  //crete and write
  {
    LxField<quad_su3> conf("conf");
    generate_hot_lx_conf(conf);
    write_ildg_gauge_conf(arg[4],conf);
  }
  
  ///////////////////////////////////////////
  
  closeNissa();
  
  return 0;
}
