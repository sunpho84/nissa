#ifndef _WRAPPERS_CU
#define _WRAPPERS_CU

#include "base/debug.hpp"

namespace cuda
{
  void crash_on_unsuccess(cudaError_t ret);
  int get_device_cont();
  int get_device();
  void set_device(int asked_dev);

  template <class T> void memcpy_to_symbol(T &dev_float,T host_float)
  {crash_on_unsuccess(cudaMemcpyToSymbol(dev_float,&host_float,sizeof(T),0,cudaMemcpyHostToDevice));}
}

#endif
