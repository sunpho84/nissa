#ifndef _WRAPPERS_CU
#define _WRAPPERS_CU

#include "base/debug.hpp"
#include "macros.hpp"

namespace cuda
{
  HOST void crash_on_unsuccess(cudaError_t ret);
  HOST int get_device_cont();
  HOST int get_device();
  HOST void set_device(int asked_dev);
  HOST void cuda_free(void *&ptr);
  template <class T> HOST void memcpy_to_symbol(T &dev_var,T &host_var)
  {crash_on_unsuccess(cudaMemcpyToSymbol(dev_var,&host_var,sizeof(T),0,cudaMemcpyHostToDevice));}
  template <class T> HOST void cuda_malloc(T *&out,size_t num)
  {crash_on_unsuccess(cudaMalloc(&out,sizeof(T)*num,0));}
}

#endif
