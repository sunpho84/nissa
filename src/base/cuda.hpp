#ifndef _CUDA_HPP
#define _CUDA_HPP

#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif
#ifndef EXTERN_CUDA

#define EXTERN_CUDA extern
 #define INIT_CUDA_TO(var)
#else
 #define INIT_CUDA_TO(var) =var
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

/// CUDA_GLOBAL is actually the cuda attribute
# define CUDA_GLOBAL __global__

#undef EXTERN_CUDA
#undef INIT_CUDA_TO

#endif
