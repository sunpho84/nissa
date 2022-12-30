#ifndef _THREADS_HPP
#define _THREADS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <cstdint>
#include <cstdio>

#include <routines/mpi_routines.hpp>

#define TO_READ(A) A=A.getReadable()
#define TO_WRITE(A) A=A.getWritable()
#define CAPTURE(A...) A

#define SERIAL_LOOP(MIN,MAX,F...)		\
  for(auto i=MIN;i<MAX;i++)			\
    F(i)

#ifdef HAVE_OPENMP

namespace nissa
{
  template <typename IMin,
	    typename IMax,
	    typename F>
  INLINE_FUNCTION
  void openmpParallelFor(const int line,
			   const char* file,
			   const IMin min,
			   const IMax max,
			   F&& f)
  {
    double initTime=0;
    extern int rank,verbosity_lv;
    const bool print=(verbosity_lv>=1// 2
		      and rank==0);
    if(print)
      {
	printf("at line %d of file %s launching openmp loop [%ld,%ld)\n",
	   line,file,(int64_t)min,(int64_t)max);
	initTime=take_time();
      }
    
#pragma omp parallel for
    for(auto i=min;i<max;i++)
      f(i);
    
    if(print)
      printf(" finished in %lg s\n",take_time()-initTime);
  }
}

# define THREADS_PARALLEL_LOOP(ARGS...) openmpParallelFor(__LINE__,__FILE__,ARGS)

#else

# define THREADS_PARALLEL_LOOP(ARGS...) SERIAL_LOOP(ARGS)

#endif

/////////////////////////////////////////////////////////////////

#ifdef USE_CUDA

namespace nissa
{
  template <typename IMin,
	    typename IMax,
	    typename F>
  __global__
  void cudaGenericKernel(const IMin min,
			 const IMax max,
			 F f) // Needs to take by value
  {
    const auto i=min+blockIdx.x*blockDim.x+threadIdx.x;
    if(i<max)
      f(i);
  }
  
  template <typename IMin,
	    typename IMax,
	    typename F>
  void cudaParallelFor(const int line,
			 const char* file,
			 const IMin min,
			 const IMax max,
		       F f) // Needs to take by value
  {
    const auto length=(max-min);
    const dim3 blockDimension(128); // to be improved
    const dim3 gridDimension((length+blockDimension.x-1)/blockDimension.x);
    
    double initTime=0;
    extern int rank,verbosity_lv;
    const bool print=(verbosity_lv>=1// 2
		      and rank==0);
    if(print)
      {
	printf("at line %d of file %s launching kernel on loop [%ld,%ld) using blocks of size %d and grid of size %d\n",
	   line,file,(int64_t)min,(int64_t)max,blockDimension.x,gridDimension.x);
	initTime=take_time();
      }
    
    if(length>0)
      {
	cudaGenericKernel<<<gridDimension,blockDimension>>>(min,max,f);
	cudaDeviceSynchronize();
      }
    
    if(print)
      printf(" finished in %lg s\n",take_time()-initTime);
  }
}

# define CUDA_PARALLEL_LOOP(ARGS...) cudaParallelFor(__LINE__,__FILE__,ARGS)

#else

# define CUDA_PARALLEL_LOOP(ARGS...) SERIAL_LOOP(ARGS)

#endif

#define HOST_PARALLEL_LOOP(MIN,MAX,CAPTURES,INDEX,A...) \
  THREADS_PARALLEL_LOOP(MIN,MAX,[CAPTURES] (const int& INDEX) mutable A)

#define DEVICE_PARALLEL_LOOP(MIN,MAX,CAPTURES,INDEX,A...) \
  CUDA_PARALLEL_LOOP(MIN,MAX,[CAPTURES] CUDA_DEVICE (const int& INDEX) mutable A)

#if THREADS_TYPE == CUDA_THREADS
# define DEFAULT_PARALLEL_LOOP DEVICE_PARALLEL_LOOP
#else
# define DEFAULT_PARALLEL_LOOP HOST_PARALLEL_LOOP
#endif
 
#define PAR(MIN,MAX,CAPTURES,A...)		\
  DEFAULT_PARALLEL_LOOP(MIN,MAX,CAPTURE(CAPTURES),A)

#endif
