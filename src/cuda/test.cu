#include <stdio.h>

#include "base/global_variables.h"

#include "cuda_wrappers.h"

namespace cuda
{
  //initialize cuda
  void init(int asked_dev)
  {
    //count device
    int ndevice=get_device_count();
    if(asked_dev<0||asked_dev>=ndevices) crash("Asked device out of [0,%d) bound",ndevices);
    
    //try to set the device and check that it has been set
    set_device(asked_dev);
    if(get_device()!=asked_dev) crash("Unsuccess setting device %d",asked_dev); 
  }
  
  void test()
  {
    init(0);
  }
}
