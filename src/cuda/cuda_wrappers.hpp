#ifndef _CUDA_WRAPPERS_CU
#define _CUDA_WRAPPERS_CU

namespace cuda
{
  int get_device();
  int get_device_count();
  void crash_on_unsuccess(cudaError_t ret);
  void set_device(int asked_dev);
}

#endif
