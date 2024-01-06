#ifndef _CUDA_HPP
#define _CUDA_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <cstdarg>
#include <routines/mpiRoutines.hpp>

namespace nissa
{
  namespace resources
  {
    inline int _nCudaDevices;
    
    inline int _iCudaDevice;
  }
  
  inline const int& nCudaDevices=resources::_nCudaDevices;
  
  inline const int& iCudaDevice=resources::_iCudaDevice;
  
#define decryptCudaError(...) \
  internalDecryptCudaError(__LINE__,__FILE__,__VA_ARGS__)
  
  inline void internalDecryptCudaError(const int& line,
				       const char *file,
				       const cudaError_t& rc,
				       const char *templ,
				       ...)
  {
    if(rc!=cudaSuccess and isMasterRank())
      {
	va_list ap;
	va_start(ap,templ);
	char mess[1024];
	vsprintf(mess,templ,ap);
	va_end(ap);
	
	internal_crash(line,file,"%s, cuda raised error: %s",mess,cudaGetErrorString(rc));
      }
  }
  
  inline void initCuda()
  {
    int nDevices;
    if(cudaGetDeviceCount(&nDevices)!=cudaSuccess)
      crash("no CUDA enabled device");
    
    printf("Number of CUDA enabled devices on rank[%ld] (%s) : %d\n",thisRank(),mpiGetProcessorName().c_str(),nDevices);
    for(int i=0;i<nDevices;i++)
      {
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp,i);
	printf(" rank %ld CUDA Enabled device %d/%d: %d.%d\n",thisRank(),i,nDevices,deviceProp.major,deviceProp.minor);
      }
    //assumes that if we are seeing multiple gpus, there are nDevices ranks to attach to each of it
    if(nDevices!=1)
      {
	resources::_iCudaDevice=thisRank()%nDevices;
	decryptCudaError(cudaSetDevice(iCudaDevice),"Unable to set device %d",iCudaDevice);
      }
  }
}

#endif
