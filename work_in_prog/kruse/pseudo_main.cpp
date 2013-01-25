#include "../../src/nissa.h"

#include "src/interface/external_interface.h"

int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  //init geometry and communication structures
  int time=8,space=4;
  init_grid(time,space);
  
  //init interface additional variables (still incomplete)
  init_additional_variables();
  
  //////////////////////////////////
  
  //read configuration and compute plaquette
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);  
  read_ildg_gauge_conf(conf,"/Users/francesco/Prace/nissa/test/data/L4T8conf");
  master_printf("Plaquette: %16.16lg\n",global_plaquette_lx_conf(conf));
  
  //initialize bgq field
  bgq_gaugefield_init();
  
  //convert
  bgq_gaugefield_transferfrom((tmlQCD_su3**)conf);
  
  //////////////////////////////////
  
  nissa_free(conf);
  
  unset_additional_variables();
  
  close_nissa();
  
  return 0;
}
