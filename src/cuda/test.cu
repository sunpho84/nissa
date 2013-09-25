#include <stdio.h>

#include "base/global_variables.h"

void cuda_set_device()
{
  //count device
  int ndevices;
  cudaError_t stat=cudaGetDeviceCount(&ndevices);
  if(stat!=cudaSuccess) printf("%s\n",cudaGetErrorString(stat));
  printf("Found: %d devices\n",ndevices);
  
  //try to set a device
  for(int idev=0;idev<ndevices;idev++)
    {
      stat=cudaSetDevice(idev); 
      if(stat!=cudaSuccess) printf("%s\n",cudaGetErrorString(stat));
      if(stat==cudaErrorInvalidDevice) 
	{ 
	  perror("cudaSetDevice returned  cudaErrorInvalidDevice"); 
	}
      int device;
      cudaGetDevice(&device); 
      printf("cudaGetDevice()=%d\n",device); 
    }
}

void cuda_test()
{
  cuda_set_device();
}
