#include "src/common.c"
#include "src/random_source_generation.c"

int main()
{
  init_test();
  
  test(test_random_source_generation(),"Random source generation");
  
  close_appretto();
  
  return 0;
}
