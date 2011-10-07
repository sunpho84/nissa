#include "src/common.c"
#include "src/plaquette_computation.c"

int main()
{
  init_test();
  
  test(test_plaquette_computation(),"Configuration loading");
  
  close_appretto();
  
  return 0;
}
