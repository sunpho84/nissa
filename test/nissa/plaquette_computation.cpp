#include "src/common.cpp"
#include "src/plaquette_computation.cpp"

int main()
{
  init_test();
  
  test(test_plaquette_computation(),"Configuration loading");
  
  close_nissa();
  
  return 0;
}
