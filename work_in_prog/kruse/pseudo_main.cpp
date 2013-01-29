#include "../../src/nissa.h"

#include "src/interface/external_interface.h"
#include "src/kruse/bgq_field.h"
#include "src/kruse/bgq_comm.h"
#include "src/kruse/bgq_spinorfield.h"
#include "src/kruse/bgq_HoppingMatrix.h"

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
  quad_su3 *eo_conf[2]={nissa_malloc("conf_e",loc_volh+bord_volh,quad_su3),nissa_malloc("conf_e",loc_volh+bord_volh,quad_su3)};
  read_ildg_gauge_conf(conf,"/Users/francesco/Prace/nissa/test/data/L4T8conf");
  split_lx_conf_into_eo_parts(eo_conf,conf);
  master_printf("Plaquette: %16.16lg\n",global_plaquette_lx_conf(conf));
  
  //initialize bgq field
  bgq_gaugefield_init();
  bgq_indices_init();
  bgq_comm_mpi_init();
  bgq_spinorfields_init();
  
  ka0=ka1=ka2=ka3=1;
  
  //convert
  bgq_gaugefield_transferfrom((tmlQCD_su3**)conf);
  
  //allocate the spinors
  spincolor *legacy_source=nissa_malloc("legacy_source",loc_volh,spincolor);
  spincolor *legacy_dest=nissa_malloc("legacy_dest",loc_volh,spincolor);
  //allocate the internal vector
  bgq_weylfield_collection *source=bgq_spinorfields_allocate(1,(spinor*)legacy_source,loc_volh);
  bgq_weylfield_collection *dest=bgq_spinorfields_allocate(1,(spinor*)legacy_dest,loc_volh);
  
  //prepare the source
  bgq_spinorfield_prepareWrite(&(source->controlblocks[0]),(tristate)EVN,ly_legacy,false);
  memset(legacy_source,0,loc_volh*sizeof(spincolor));
  legacy_source[0][0][0][0]=1;
  
  int ip=loceo_neighup[0][0][0];
  
  //test
  spincolor *test_dest=nissa_malloc("legacy_dest",loc_volh,spincolor);
  tmn2Doe_eos(test_dest,eo_conf,legacy_source);
  master_printf("%lg\n",test_dest[ip][0][0][0]);
  
  bgq_HoppingMatrix(ODD,&(dest->controlblocks[0]),&(source->controlblocks[0]),hm_nokamul);
  bgq_spinorfield_prepareRead(&(dest->controlblocks[0]),(tristate)ODD,false,false,false,false,true);
  master_printf("%lg\n",legacy_dest[ip][0][0][0]);
  
  //////////////////////////////////
  
  nissa_free(conf);
  nissa_free(legacy_source);
  
  unset_additional_variables();
  
  close_nissa();
  
  return 0;
}
