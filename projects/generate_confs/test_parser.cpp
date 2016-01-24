#include "nissa.hpp"

#include "driver.hpp"

using namespace nissa;

void in_main(int narg,char **arg)
{
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //open input file
  open_input(arg[1]);
  
  driver_t *driver=new driver_t(input_global);
  parser_parse(driver);
  
  driver->master_fprintf(stdout);
  
  master_printf("read %d quarks\n",(int)driver->quarks.size());
  master_printf("nlevels %d\n",(int)driver->stout_pars.nlevels);
  
  delete driver;
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}

