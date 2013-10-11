#include <stdio.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "routines/ios.hpp"

#include "fields.hpp"
#include "global_variables.hpp"
#include "new_types.hpp"
#include "wrappers.hpp"

namespace cuda
{
  //check that we can start cuda
  void check_startable()
  {
    //check that grid and e/o geometry are inited
    if(!nissa::grid_inited) nissa::crash("Nissa grid not inited");
    if(!nissa::eo_geom_inited) nissa::crash("e/o geom not inited");
  }
  
  //initialize cuda
  void init(int asked_dev)
  {
    //perform some check and select the device
    check_startable();
    set_device(asked_dev);
    
    //copy geometry information
    nissa::master_printf("loc_vol: %d\n",nissa::loc_vol);
    memcpy_to_symbol(loc_vol,nissa::loc_vol);
    nissa::master_printf("loc_volh: %d\n",nissa::loc_volh);
    //memcpy_to_symbol(loc_volh,nissa::loc_volh);
    
    //define the field
    float_gauge_field test;
  }
  
  void test()
  {
    init(1);
  }
}
