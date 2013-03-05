#include "nissa.h"

void in_main(int narg,char **arg)
{
  init_grid(16,16);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  
  close_nissa();
  
  return 0;
}
