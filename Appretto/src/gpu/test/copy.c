#include <stdio.h>

#include "gpu.h"

int rank=0;

#define size 14

int main()
{
  char cpu1[size]={"Ciao mondo!\n"};
  char *gpu1=gpu_allocate_vector(size,"ciao_mondo1");
  char *gpu2=gpu_allocate_vector(size,"ciao_mondo2");
  char cpu2[size]={""};


  memcpy_cpu_to_gpu(gpu1,cpu1,size);
  memcpy_gpu_to_gpu(gpu2,gpu1,size);
  memcpy_gpu_to_cpu(cpu2,gpu2,size);
  
  printf("%s",cpu2);
  
  gpu_free(gpu1);
  gpu_free(gpu2);
  
  return 0;
}
