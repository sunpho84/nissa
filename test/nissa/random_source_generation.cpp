#include "src/common.cpp"
#include "src/random_source_generation.cpp"

int main()
{
  init_test();
  
  test(test_random_source_generation(),"Random source generation");
  
  close_nissa();
  
  return 0;
}
