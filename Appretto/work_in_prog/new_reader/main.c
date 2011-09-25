#include "appretto.h"


int main(int narg,char **arg)
{
  init_appretto();
  if(narg<2) crash("Use: %s input_file",arg[0]);

  open_input(arg[1]);

  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));

  char filename[1024];
  read_str_str("Filename",filename,1024);

  close_input();
  
  init_grid();
  
  quad_su3 *conf=appretto_malloc("conf",loc_vol,quad_su3);
  
  read_gauge_conf(conf,filename);
  
  appretto_free(conf);
  
  close_appretto();
  
  return 0;
}
