#ifndef _CUDA_HPP
#define _CUDA_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/memory_manager.hpp>
#include <routines/mpi_routines.hpp>

#include <geometry/geometry_lx.hpp>

namespace nissa
{
  /// Number of cuda device
  inline int nCudaDevices;
  
  /// Id of the current cuda device
  inline int iCudaDevice;
  
  /// Properties of the curren device
  inline cudaDeviceProp deviceProp;
  
  /// Initialize cuda
  inline void init_cuda()
  {
    if(cudaGetDeviceCount(&nCudaDevices)!=cudaSuccess)
      CRASH("no CUDA enabled device");
    
    printf("Number of CUDA enabled devices on rank[%d] (%s) : %d\n",rank,MPI_get_processor_name().c_str(),nCudaDevices);
    for(int i=0;i<nCudaDevices;i++)
      {
	cudaGetDeviceProperties(&deviceProp,i);
	printf(" rank %d CUDA Enabled device %d/%d: %d.%d\n",rank,i,nCudaDevices,deviceProp.major,deviceProp.minor);
      }
    //assumes that if we are seeing multiple gpus, there are nDevices ranks to attach to each of it
    if(nCudaDevices!=1)
      {
	iCudaDevice=rank%nCudaDevices;
	DECRYPT_CUDA_ERROR(cudaSetDevice(iCudaDevice),"Unable to set device %d",iCudaDevice);
      }
  }
  
  /// Returns the cudaMemcpyKind corresponding to the required copy from/to
  template <MemorySpace TO,
	    MemorySpace FROM>
  inline cudaMemcpyKind memcpyKindForCopy;
  
  /// Returns the cudaMemcpyKind corresponding to the CPU-CPU copy
  template <>
  inline cudaMemcpyKind memcpyKindForCopy<MemorySpace::CPU,MemorySpace::CPU> =
    cudaMemcpyHostToHost;
  
  /// Returns the cudaMemcpyKind corresponding to the CPU-GPU copy
  template <>
  inline cudaMemcpyKind memcpyKindForCopy<MemorySpace::GPU,MemorySpace::CPU> =
    cudaMemcpyHostToDevice;
  
  /// Returns the cudaMemcpyKind corresponding to the GPU-CPU copy
  template <>
  inline cudaMemcpyKind memcpyKindForCopy<MemorySpace::CPU,MemorySpace::GPU> =
    cudaMemcpyDeviceToHost;
  
  /// Returns the cudaMemcpyKind corresponding to the GPU-GPU copy
  template <>
  inline cudaMemcpyKind memcpyKindForCopy<MemorySpace::GPU,MemorySpace::GPU> =
    cudaMemcpyDeviceToDevice;
}

#endif
