#include <stdio.h>

#include "base/global_variables.hpp"
#include "base/debug.hpp"

#include "cuda_wrappers.hpp"

namespace cuda
{
  //initialize cuda
  void init(int asked_dev)
  {
    //count device
    int ndevices=get_device_count();
    if(asked_dev<0||asked_dev>=ndevices) nissa::crash("Asked device out of [0,%d) bound",ndevices);
    
    //try to set the device and check that it has been set
    set_device(asked_dev);
    if(get_device()!=asked_dev) nissa::crash("Unsuccess setting device %d",asked_dev); 
  }
  
  void test()
  {
    init(0);
  }
}
