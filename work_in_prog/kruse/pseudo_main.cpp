#include "../../src/nissa.h"

#include "src/interface/external_interface.h"

int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  init_grid(8,4);
  
  init_additional_variables();
  
  //////////////////////////////////
  
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  
  read_ildg_gauge_conf(conf,"/Users/francesco/Prace/nissa/test/data/L4T8conf");
  master_printf("Plaquette: %16.16lg\n",global_plaquette_lx_conf(conf));
  
  bgq_gaugefield_init();
  
  nissa_free(conf);
  
  //////////////////////////////////
  
  unset_additional_variables();
  
  close_nissa();
  
  return 0;
}
