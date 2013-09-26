#include "base/debug.hpp"

namespace cuda
{
  //crash translating cuda error
  void crash_on_unsuccess(cudaError_t ret)
  {if(ret!=cudaSuccess) nissa::crash("%s",cudaGetErrorString(ret));}
  
  //return the number of devices
  int get_device_count()
  {
    int ndevices;
    crash_on_unsuccess(cudaGetDeviceCount(&ndevices));
    
    return ndevices;
  }
  
  //return used devices
  int get_device()
  {
    int idev;
    crash_on_unsuccess(cudaGetDevice(&idev));
    
    return idev;
  }
  
  //set the device to be used
  void set_device(int asked_dev)
  {
    //count devices and check that we are asking for something present
    int ndevices=get_device_count();
    if(asked_dev<0||asked_dev>=ndevices) nissa::crash("Asked device out of [0,%d) bound",ndevices);
    
    //activate and check
    crash_on_unsuccess(cudaSetDevice(asked_dev));
    if(get_device()!=asked_dev) nissa::crash("Unsuccess setting device %d",asked_dev);
  }
}
