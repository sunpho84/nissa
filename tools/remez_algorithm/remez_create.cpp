#include "nissa.hpp"

using namespace nissa;

int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  rat_approx_t rat;
  generate_approx_of_maxerr(rat,pow(0.0027,2),10,1e-3,-1,32);
  
  close_nissa();
  
  return 0;
}
