#ifndef _CUDA_HPP
#define _CUDA_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif
#ifndef EXTERN_CUDA

#define EXTERN_CUDA extern
 #define INIT_TO(var)
#else
 #define INIT_TO(var) =var
#endif

namespace nissa
{
#ifdef USE_CUDA
  EXTERN_CUDA int nCudaDevices;
  EXTERN_CUDA int iCudaDevice;
#endif
  
  void init_cuda();
  
  /// True or false depending on whether we are compiling on device
  [[ maybe_unused ]]
  constexpr bool CompilingForDevice=
#ifdef COMPILING_FOR_DEVICE
    true
#else
    false
#endif
    ;
}

#ifdef USE_CUDA

 /// CUDA_DEVICE is actually the cuda attribute
# define CUDA_DEVICE __device__
 
 /// CUDA_GLOBAL is actually the cuda attribute
# define CUDA_GLOBAL __global__
 
 /// CUDA_HOST is actually the cuda attribute
# define CUDA_HOST __host__
 
#else
 
 /// CUDA_DEVICE is a dummy macro
# define CUDA_DEVICE
 
 /// CUDA_HOST is a dummy macro
# define CUDA_HOST
 
 /// CUDA_GLOBAL is a dummy macro
# define CUDA_GLOBAL
 
#endif

/// Put together CUDA_HOST and CUDA_DEVICE
#define CUDA_HOST_DEVICE CUDA_HOST CUDA_DEVICE

#undef EXTERN_CUDA
#undef INIT_TO

#endif
