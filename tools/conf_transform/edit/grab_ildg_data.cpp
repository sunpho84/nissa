#include "nissa.hpp"

using namespace nissa;

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<5) crash("use: %s L T file_in file_out",arg[0]);
  
  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  
  //Init the MPI grid
  init_grid(T,L);
  //////////////////////////// read the conf /////////////////////////////
  
  quad_su3 *conf=nissa_malloc("conf",locVolWithBord.nastyConvert(),quad_su3);
  read_ildg_gauge_conf(conf,arg[3]);
  
  //////////////////////////// write the conf ////////////////////////////
  
  write_ildg_gauge_conf(arg[4],conf,64);
  nissa_free(conf);
  
  ///////////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
