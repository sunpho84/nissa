#include <cuda.h>
#include <stdio.h>
#include <cuda_runtime.h>

extern "C" void choose_gpu()
{
  //check that we have at least one gpu
  if(find_gpu()==0)
    {
      fprintf(stderr,"Error: no gpu found.\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  
  if(cudaSetDevice(0)!=cudaSuccess)
    {
      fprintf(stderr,"Could not set active device.\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }  
  
  int dev_num;
  cudaGetDevice(&dev_num);
  printf("Rank %d choose device: %d\n",rank,dev_num);
}
