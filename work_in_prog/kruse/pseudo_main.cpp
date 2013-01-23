#include "src/interface/interface.h"

#include "../../src/nissa.h"

int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  init_grid(4,4);
  
  init_additional_variables();
  
  unset_additional_variables();
  
  close_nissa();
  
  return 0;
}
