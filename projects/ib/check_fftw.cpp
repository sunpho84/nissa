#include <nissa.hpp>

using namespace nissa;

void in_main(int narg,char **arg)
{
    
  glb_size[0]=glb_size[1]=glb_size[2]=2;
  glb_size[3]=6;
  init_grid(-1,-1);
  
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}

