#ifndef _CUDA_MACROS_HPP
#define _CUDA_MACROS_HPP

#define GLOBAL __global__
#define DEVICE __device__
#define HOST __host__
#define CONSTANT __constant__
#define DEVICE_CONSTANT DEVICE CONSTANT

#define NCUDA_THREADS 128

#endif
