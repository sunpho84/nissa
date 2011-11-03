#include "src/common.c"
#include "src/write_and_read.c"

int main()
{
  init_test();
  
  test(test_write_and_read(),"Write and read");
  
  close_appretto();
  
  return 0;
}
