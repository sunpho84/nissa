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
  {crash_on_unsuccess(cudaSetDevice(asked_dev));}
  
  //copy a float to the symbols
  void memcpy_to_symbol(float dev_float,float host_float)
  {crash_on_unsuccess(cudaMemcpyToSymbol(dev_float,&host_float,sizeof(float),0,cudaMemcpyHostToDevice));}
}
