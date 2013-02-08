#include "../../src/nissa.h"

#include "src/interface/external_interface.h"
#include "src/kruse/bgq_field.h"
#include "src/kruse/bgq_comm.h"
#include "src/kruse/bgq_spinorfield.h"
#include "src/kruse/bgq_HoppingMatrix.h"

int bgq_field_initted=false;
bgq_weylfield_collection *bgq_source,*bgq_temp,*bgq_dest;

//allocate the internal vector
void allocate_bgq_field(spincolor *dest,spincolor *temp,spincolor *source)
{
  bgq_field_initted=true;
  
  bgq_source=bgq_spinorfields_allocate(1,(spinor*)source,loc_volh);
  bgq_temp=bgq_spinorfields_allocate(1,(spinor*)temp,loc_volh);
  bgq_dest=bgq_spinorfields_allocate(1,(spinor*)dest,loc_volh);
}

void app(spincolor *legacy_dest,spincolor *legacy_temp,spincolor *legacy_source,spincolor *ext_source=NULL)
{
  if(!bgq_field_initted) allocate_bgq_field(legacy_dest,legacy_temp,legacy_source);

  //if an ext source is passed, copy it
  if(ext_source!=NULL)
    {
      //prepare the source to be wrote
      bgq_spinorfield_prepareWrite(&(bgq_source->controlblocks[0]),(tristate)EVN,ly_legacy,false);
      //write inside
      vector_copy(legacy_source,ext_source);
    }
  
  //apply OE
  bgq_HoppingMatrix(ODD,&(bgq_temp->controlblocks[0]),&(bgq_source->controlblocks[0]),hm_nokamul);

  //apply EO
  bgq_HoppingMatrix(EVN,&(bgq_dest->controlblocks[0]),&(bgq_temp->controlblocks[0]),hm_nokamul);
}

void init(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  //init geometry and communication structures
  int time=8,space=4;
  init_grid(time,space);
  
  //init interface
  init_interface();
}

int main(int narg,char **arg)
{
  init(narg,arg);
  
  //////////////////////////////////
  
  //read configuration and compute plaquette
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);  
  quad_su3 *eo_conf[2]={nissa_malloc("conf_e",loc_volh+bord_volh,quad_su3),nissa_malloc("conf_e",loc_volh+bord_volh,quad_su3)};
  read_ildg_gauge_conf(conf,"/Users/francesco/Prace/nissa/test/data/L4T8conf");
  //vector_reset(conf);
  //nissa_loc_vol_loop(ivol)
  //for(int i=0;i<4;i++)
  //su3_put_to_id(conf[ivol][i]);
  
  split_lx_conf_into_eo_parts(eo_conf,conf);
  master_printf("Plaquette: %16.16lg\n",global_plaquette_lx_conf(conf));
  
  //convert
  bgq_gaugefield_transferfrom((su3*)conf);
  
  //allocate the spinors
  spincolor *legacy_source=nissa_malloc("legacy_source",loc_volh+bord_volh,spincolor);
  spincolor *legacy_temp=nissa_malloc("legacy_temp",loc_volh+bord_volh,spincolor);
  spincolor *legacy_dest=nissa_malloc("legacy_dest",loc_volh,spincolor);
  spincolor *ori_source=nissa_malloc("ori_source",loc_volh+bord_volh,spincolor);
  
  //prepare the source
  vector_reset(ori_source);
  ori_source[0][0][0][0]=1;
  
  //vector_reset(ori_source);
  //nissa_loc_volh_loop(ivol)
  //for(int mu=0;mu<4;mu++)
  //ori_source[ivol][mu/3][mu%3][0]=glb_coord_of_loclx[loclx_of_loceo[EVN][ivol]][mu];
  
  //test
  tmn2Doe_eos(legacy_temp,eo_conf,ori_source);
  tmn2Deo_eos(legacy_dest,eo_conf,legacy_temp);
  
  for(int mu=0;mu<4;mu++)
    {
      int ip=loceo_neighup[ODD][0][mu];
      master_printf("leg %d ip %d %lg\n",mu,ip,legacy_dest[ip][0][0][0]);
      ip=loceo_neighdw[ODD][0][mu];
      master_printf("leg %d, ip %d %lg\n",mu,ip,legacy_dest[ip][0][0][0]);
    }
  
  app(legacy_dest,legacy_temp,legacy_source,ori_source);
  bgq_spinorfield_prepareRead(&(bgq_temp->controlblocks[0]),(tristate)ODD,false,false,false,false,true);
  bgq_spinorfield_prepareRead(&(bgq_dest->controlblocks[0]),(tristate)EVN,false,false,false,false,true);
  
  for(int mu=0;mu<4;mu++)
    {
      int ip=loceo_neighup[ODD][0][mu];
      master_printf("new %d ip %d %lg\n",mu,ip,legacy_dest[ip][0][0][0]);
      ip=loceo_neighdw[ODD][0][mu];
      master_printf("new %d ip %d %lg\n",mu,ip,legacy_dest[ip][0][0][0]);
    }
  
  //////////////////////////////////
  
  nissa_free(conf);
  nissa_free(legacy_source);
  
  unset_interface();
  
  close_nissa();
  
  return 0;
}
