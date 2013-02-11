#include "nissa.h"

#include "spi.h"

int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  init_grid(16,16);
  
  init_SPI();
  
  test_SPI_comm();
  
  close_nissa();
  
  return 0;
}
