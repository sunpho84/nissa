#ifndef _CUDA_HPP
#define _CUDA_HPP

#ifndef EXTERN_CUDA
 #define EXTERN_CUDA extern
 #define INIT_TO(var)
#else
 #define INIT_TO(var) =var
#endif

#include <cuda_runtime.h>

namespace nissa
{
  EXTERN_CUDA int nCudaDevices;
  EXTERN_CUDA int iCudaDevice;
  
  void init_cuda();
}

#undef EXTERN_CUDA
#undef INIT_TO

#endif
