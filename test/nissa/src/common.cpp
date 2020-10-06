#pragma once

#include "nissa.h"

void init_test()
{
  init_nissa();
  
  //init the grid
  init_grid(8,4);
}

void test(int passed,const char *test_name)
{
  if(!passed) master_printf("\n%s test not passed!\n\n",test_name);
  else master_printf("\n%s test passed\n\n",test_name);
  
  master_printf("################ %s test finished ###############\n",test_name);
}

