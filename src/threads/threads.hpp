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
  
#ifdef USE_OPENMP
  
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
  
#endif
  
  /////////////////////////////////////////////////////////////////
  
#ifdef USE_CUDA
  
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
	decript_cuda_error(cudaDeviceSynchronize(),"during kernel execution");
      }
    
    if(print)
      printf(" finished in %lg s\n",take_time()-initTime);
  }
#endif
  
  // End of backend for parallel for
  
  /////////////////////////////////////////////////////////////////
  
  // Associate host or device loop with the appropriate backend
  
  // At the moment single backend is availbale for both sides, so no
  // choice is possible, but in future one might want to select the
  // backend from the configuration command
  
#ifdef USE_OPENMP
# define HOST_PARALLEL_LOOP_ROUTINE openmpParallelFor
#else
# define HOST_PARALLEL_LOOP_ROUTINE serialFor
#endif
  
  /////////////////////////////////////////////////////////////////
  
#ifdef USE_CUDA
# define DEVICE_PARALLEL_LOOP_ROUTINE cudaParallelFor
#else
# define DEVICE_PARALLEL_LOOP_ROUTINE HOST_PARALLEL_LOOP_ROUTINE
#endif

#define HOST_PARALLEL_LOOP(MIN,MAX,CAPTURES,INDEX,BODY...)		\
  HOST_PARALLEL_LOOP_ROUTINE(__LINE__,__FILE__,MIN,MAX,[CAPTURES] (const auto& INDEX) MUTABLE_INLINE_ATTRIBUTE BODY)

#define DEVICE_PARALLEL_LOOP(MIN,MAX,CAPTURES,INDEX,BODY...)		\
  DEVICE_PARALLEL_LOOP_ROUTINE(__LINE__,__FILE__,MIN,MAX,[CAPTURES] CUDA_DEVICE (const auto& INDEX) mutable BODY)
  
  ///////////////////////////////////////////////////////////////// this part must go
  
#if USE_CUDA
# define DEFAULT_PARALLEL_LOOP DEVICE_PARALLEL_LOOP
#else
# define DEFAULT_PARALLEL_LOOP HOST_PARALLEL_LOOP
#endif
 
#define PAR(MIN,MAX,CAPTURES,BODY...)			\
  DEFAULT_PARALLEL_LOOP(MIN,MAX,CAPTURE(CAPTURES),BODY)
}

#endif
