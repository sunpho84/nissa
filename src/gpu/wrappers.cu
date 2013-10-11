#include "base/debug.hpp"
#include "routines/ios.hpp"
#include "macros.hpp"

namespace cuda
{
  //crash translating cuda error
  HOST void crash_on_unsuccess(cudaError_t ret)
  {if(ret!=cudaSuccess) nissa::crash("%s",cudaGetErrorString(ret));}
  
  //return the number of devices
  HOST int get_device_count()
  {
    int ndevices;
    crash_on_unsuccess(cudaGetDeviceCount(&ndevices));
    
    return ndevices;
  }
  
  //free on the device
  HOST void cuda_free(void *&ptr)
  {
    crash_on_unsuccess(cudaFree(ptr));
    ptr=NULL;
  }
  
  //return used devices
  HOST int get_device()
  {
    int idev;
    crash_on_unsuccess(cudaGetDevice(&idev));
    
    return idev;
  }
  
  //set the device to be used
  HOST void set_device(int asked_dev)
  {
    //count devices and check that we are asking for something present
    int ndevices=get_device_count();
    nissa::master_printf("Number of present devices: %d\n",ndevices);
    if(asked_dev<0||asked_dev>=ndevices) nissa::crash("Asked device out of [0,%d) bound",ndevices);
    
    //activate and check
    crash_on_unsuccess(cudaSetDevice(asked_dev));
    if(get_device()!=asked_dev) nissa::crash("Unsuccess setting device %d",asked_dev);
    nissa::master_printf("Correctly set device: %d\n",asked_dev);
  }
}
