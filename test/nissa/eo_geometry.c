#include "src/common.c"
#include "src/eo_geometry.c"

int main()
{
  init_test();
  
  test(test_eo_geometry(),"E/O geometry");
  
  close_nissa();
  
  return 0;
}
