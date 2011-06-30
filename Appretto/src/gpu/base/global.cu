#pragma once

//duplicated from cpu
char *cpu_allocate_vector(int length,char *tag)
{
  char *out=(char*)malloc(length);
  if(out==NULL && rank==0)
    {
      fprintf(stderr,"Error during allocation of %s\n",tag);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  return out;
}
su3c *cpu_allocate_su3c(int length,char *tag){return (su3c*)cpu_allocate_vector(length*sizeof(su3c),tag);}
su3f *cpu_allocate_su3f(int length,char *tag){return (su3f*)cpu_allocate_vector(length*sizeof(su3f),tag);}
quad_su3c *cpu_allocate_quad_su3c(int length,char *tag){return (quad_su3c*)cpu_allocate_vector(length*sizeof(quad_su3c),tag);}
quad_su3f *cpu_allocate_quad_su3f(int length,char *tag){return (quad_su3f*)cpu_allocate_vector(length*sizeof(quad_su3f),tag);}

//allocate vectors of the required length
char *gpu_allocate_vector(int length,char *tag)
{
  char *out;
  
  if(cudaMalloc((void**)&out,length)!=cudaSuccess && rank==0)
    {
      fprintf(stderr,"Error during allocation of %s\n",tag);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  
  return out;
}
su3c *gpu_allocate_su3c(int length,char *tag){return (su3c*)gpu_allocate_vector(length*sizeof(su3c),tag);}
su3f *gpu_allocate_su3f(int length,char *tag){return (su3f*)gpu_allocate_vector(length*sizeof(su3f),tag);}
quad_su3c *gpu_allocate_quad_su3c(int length,char *tag){return (quad_su3c*)gpu_allocate_vector(length*sizeof(quad_su3c),tag);}
quad_su3f *gpu_allocate_quad_su3f(int length,char *tag){return (quad_su3f*)gpu_allocate_vector(length*sizeof(quad_su3f),tag);}

//free a vector on the gpu
extern "C" void gpu_free(void *out,char *tag)
{
  if(out!=NULL) cudaFree(out);
  else
    if(rank==0)
      {
	fprintf(stderr,"Error during de-allocation of %s, already free\n",tag);
	MPI_Abort(MPI_COMM_WORLD,1);
      }
}

void memcpy_cgpu_to_cgpu(void *dest,const void *source,int size,cudaMemcpyKind dir)
{
  if(cudaMemcpy(dest,source,size,dir)!=cudaSuccess && rank==0)
    {
      fprintf(stderr,"Error during memory copy!\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
}

extern "C" void memcpy_cpu_to_gpu(void *gpu,const void *cpu,int size)
{memcpy_cgpu_to_cgpu(gpu,cpu,size,cudaMemcpyHostToDevice);}

extern "C" void memcpy_gpu_to_cpu(void *cpu,const void *gpu,int size)
{memcpy_cgpu_to_cgpu(cpu,gpu,size,cudaMemcpyDeviceToHost);}

extern "C" void memcpy_gpu_to_gpu(void *gpu_dest,const void *gpu_source,int size)
{memcpy_cgpu_to_cgpu(gpu_dest,gpu_source,size,cudaMemcpyDeviceToDevice);}

