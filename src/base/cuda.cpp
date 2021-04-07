#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_CUDA
 #include <base/cuda.hpp>

#include "routines/mpi_routines.hpp"

#include "geometry/geometry_lx.hpp"

namespace nissa
{
  void init_cuda()
  {
    int nDevices;
    if(cudaGetDeviceCount(&nDevices)!=cudaSuccess)
      crash("no CUDA enabled device");
    
    printf("Number of CUDA enabled devices on rank[%d] (%s) : %d\n",rank,MPI_get_processor_name().c_str(),nDevices);
    for(int i=0;i<nDevices;i++)
      {
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp,i);
	printf(" rank %d CUDA Enabled device %d/%d: %d.%d\n",rank,i,nDevices,deviceProp.major,deviceProp.minor);
      }
    
    //assumes that if we are seeing multiple gpus, there are nDevices ranks to attach to each of it
    if(nDevices!=1)
      {
	iCudaDevice=rank%nDevices;
	decript_cuda_error(cudaSetDevice(iCudaDevice),"Unable to set device %d",iCudaDevice);
      }
  }
}
