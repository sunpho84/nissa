#include "src/common.cpp"
#include "src/write_and_read.cpp"

int main()
{
  init_test();
  
  test(test_write_and_read(),"Write and read");
  
  close_nissa();
  
  return 0;
}
