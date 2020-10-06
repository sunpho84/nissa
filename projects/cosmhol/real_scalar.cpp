#include "nissa.hpp"

using namespace nissa;

void in_main(int narg,char **arg)
{
  if(narg<2) crash("use %s nx",arg[0]);
  int X=atoi(arg[1]);
  int T=2*X,L=X;
  init_grid(T,L);
  glb_size[0]=T;
  for(int mu=1;mu<4;mu++) glb_size[mu]=L;
  
  start_loc_rnd_gen(1000);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
