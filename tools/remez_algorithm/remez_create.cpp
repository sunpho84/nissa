#include "nissa.hpp"

using namespace nissa;

int main(int narg,char **arg)
{
  initNissa(narg,arg);
  
  rat_approx_t rat;
  generate_approx_of_maxerr(rat,2.3273336186266341e-07,1,7.418800292288640e-09,1,8);
  
  rat.master_fprintf(stdout);
  closeNissa();
  
  return 0;
}
