#ifndef _CUDA_HPP
#define _CUDA_HPP

#ifndef EXTERN_CUDA
 #define EXTERN_CUDA extern
 #define INIT_TO(var)
#else
 #define INIT_TO(var) =var
#endif

namespace nissa
{
  EXTERN_CUDA int nCudaDevices;
  EXTERN_CUDA int iCudaDevice;
  
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

#undef EXTERN_CUDA
#undef INIT_TO

#endif
