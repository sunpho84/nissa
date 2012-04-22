#include "src/common.cpp"
#include "src/eo_geometry.cpp"

int main()
{
  init_test();
  
  test(test_eo_geometry(),"E/O geometry");
  
  close_nissa();
  
  return 0;
}
