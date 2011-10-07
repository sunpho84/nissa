#pragma once

#include "appretto.h"

void init_test()
{
  init_appretto();
  
  //set the lattice grid
  glb_size[0]=8;
  glb_size[1]=4;
  
  //init the grid
  init_grid();
}

void test(int passed,const char *test_name)
{
  if(!passed) master_printf("\n%s test not passed!\n\n",test_name);
  else master_printf("\n%s test passed\n\n",test_name);
}

