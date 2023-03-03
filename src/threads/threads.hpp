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

#if THREADS_TYPE == OPENMP_THREADS

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
	printf("at line %d of file %s launching openmp loop [%d,%d)\n",
	   line,file,(int)min,(int)max);
	initTime=take_time();
      }
    
#pragma omp parallel for
    for(int i=min;i<(int)max;i++)
      f(IMax(i));
    
    if(print)
      printf(" finished in %lg s\n",take_time()-initTime);
  }
}

# define THREADS_PARALLEL_LOOP(ARGS...) openmpParallelFor(__LINE__,__FILE__,ARGS)

#else

namespace nissa
{
  template <typename IMin,
        typename IMax,
        typename F>
  INLINE_FUNCTION
  void serialFor(const int line,
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
    printf("at line %d of file %s launching serial loop [%d,%d)\n",
       line,file,(int)min,(int)max);
    initTime=take_time();
      }

    for(int i=min;i<(int)max;i++)
      f(IMax(i));

    if(print)
      printf(" finished in %lg s\n",take_time()-initTime);
  }
}

# define THREADS_PARALLEL_LOOP(ARGS...) serialFor(__LINE__,__FILE__,ARGS)

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
    const IMax i=(int)(min+blockIdx.x*blockDim.x+threadIdx.x);
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
    const auto length=(int)max-(int)min;
    const dim3 blockDimension(128); // to be improved
    const dim3 gridDimension((length+blockDimension.x-1)/blockDimension.x);
    
    double initTime=0;
    extern int rank,verbosity_lv;
    const bool print=(verbosity_lv>=1// 2
		      and rank==0);
    if(print)
      {
	printf("at line %d of file %s launching kernel on loop [%d,%d) using blocks of size %d and grid of size %d\n",
	   line,file,(int)min,(int)max,blockDimension.x,gridDimension.x);
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

# define CUDA_PARALLEL_LOOP(ARGS...) serialFor(__LINE__,__FILE__,ARGS)

#endif

#define HOST_PARALLEL_LOOP(MIN,MAX,CAPTURES,INDEX,A...) \
  THREADS_PARALLEL_LOOP(MIN,MAX,[CAPTURES] (const auto& INDEX) MUTABLE_INLINE_ATTRIBUTE A)

#define DEVICE_PARALLEL_LOOP(MIN,MAX,CAPTURES,INDEX,A...) \
  CUDA_PARALLEL_LOOP(MIN,MAX,[CAPTURES] CUDA_DEVICE (const auto& INDEX) mutable A)

#if THREADS_TYPE == CUDA_THREADS
# define DEFAULT_PARALLEL_LOOP DEVICE_PARALLEL_LOOP
#else
# define DEFAULT_PARALLEL_LOOP HOST_PARALLEL_LOOP
#endif
 
#define PAR(MIN,MAX,CAPTURES,A...)		\
  DEFAULT_PARALLEL_LOOP(MIN,MAX,CAPTURE(CAPTURES),A)

#endif
