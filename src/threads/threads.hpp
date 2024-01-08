#ifndef _THREADS_HPP
#define _THREADS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <cstdint>
#include <cstdio>

#ifdef ENABLE_DEVICE_CODE
# include <base/cuda.hpp>
#endif
#include <expr/comp.hpp>
#include <routines/mpiRoutines.hpp>

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
    using Idx=std::common_type_t<IMin,IMax>;
    
    double initTime=0;
    extern int verbosity_lv;
    const bool print=(verbosity_lv>=1// 2
		      and isMasterRank());
    if(print)
      {
	printf("at line %d of file %s launching serial loop [%d,%d)\n",
	       line,file,(int)min,(int)max);
	initTime=takeTime();
      }
    
    for(auto i=compDecay(min);i<compDecay(max);i++)
      f(Idx(i));
    
    if(print)
      printf(" finished in %lg s\n",takeTime()-initTime);
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
    using Idx=
      std::common_type_t<IMin,IMax>;
    
    double initTime=0;
    extern int verbosityLv;
    const bool print=(verbosityLv>=1// 2
		      and isMasterRank());
    if(print)
      {
	printf("at line %d of file %s launching openmp loop [%d,%d)\n",
	       line,file,(int)compDecay(min),(int)compDecay(max));
	initTime=takeTime();
      }
    
#pragma omp parallel for
    for(std::decay_t<decltype(compDecay(std::declval<Idx>()))> i=compDecay(min);
	i<compDecay(max);
	i++)
      f(Idx(i));
    
    if(print)
      printf(" finished in %lg s\n",takeTime()-initTime);
  }
  
#endif
  
  /////////////////////////////////////////////////////////////////
  
#ifdef ENABLE_DEVICE_CODE
  
  template <typename IMin,
	    typename IMax,
	    typename F>
  __global__
  void cudaGenericKernel(const IMin min,
			 const IMax max,
			 F f) // Needs to take by value
  {
    using Idx=
      std::common_type_t<IMin,IMax>;
    
    const Idx i=
      min+blockIdx.x*blockDim.x+threadIdx.x;
    
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
		       F&& f)
  {
    const auto length=compDecay(max)-compDecay(min);
    const dim3 blockDimension(128); // to be improved
    const dim3 gridDimension((length+blockDimension.x-1)/blockDimension.x);
    
    double initTime=0;
    extern int verbosity_lv;
    const bool print=(verbosity_lv>=1// 2
		      and isMasterRank());
    if(print)
      {
	printf("at line %d of file %s launching kernel on loop [%d,%d) using blocks of size %d and grid of size %d\n",
	       line,file,(int)compDecay(min),(int)compDecay(max),blockDimension.x,gridDimension.x);
	initTime=take_time();
      }
    
    if(length>0)
      {
	cudaGenericKernel<<<gridDimension,blockDimension>>>(min,max,std::forward<F>(f));
	decryptCudaError(cudaDeviceSynchronize(),"during kernel execution");
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
  
#ifdef ENABLE_DEVICE_CODE
# define DEVICE_PARALLEL_LOOP_ROUTINE cudaParallelFor
#else
# define DEVICE_PARALLEL_LOOP_ROUTINE HOST_PARALLEL_LOOP_ROUTINE
#endif

#define HOST_PARALLEL_LOOP(MIN,MAX,CAPTURES,INDEX,BODY...)		\
  HOST_PARALLEL_LOOP_ROUTINE(__LINE__,__FILE__,MIN,MAX,[CAPTURES] (const auto& INDEX) MUTABLE_INLINE_ATTRIBUTE BODY)

#define DEVICE_PARALLEL_LOOP(MIN,MAX,CAPTURES,INDEX,BODY...)		\
  DEVICE_PARALLEL_LOOP_ROUTINE(__LINE__,__FILE__,MIN,MAX,[CAPTURES] DEVICE_ATTRIB (const auto& INDEX) mutable BODY)
  
  ///////////////////////////////////////////////////////////////// this part must go
  
#ifdef ENABLE_DEVICE_CODE
# define DEFAULT_PARALLEL_LOOP DEVICE_PARALLEL_LOOP
#else
# define DEFAULT_PARALLEL_LOOP HOST_PARALLEL_LOOP
#endif
 
#define PAR(MIN,MAX,CAPTURES,INDEX,BODY...)		\
  DEFAULT_PARALLEL_LOOP(MIN,MAX,CAPTURE(CAPTURES),INDEX,BODY)

#ifdef ENABLE_DEVICE_CODE
# define PAR_ON_EXEC_SPACE(EXEC_SPACE,MIN,MAX,CAPTURES,INDEX,BODY...)	\
  static_assert(EXEC_SPACE.hasUniqueExecSpace(),"Not uinique exec space!"); \
  if constexpr(EXEC_SPACE==execOnCPU)					\
    HOST_PARALLEL_LOOP(MIN,MAX,CAPTURE(CAPTURES),INDEX,BODY);		\
  else									\
    DEVICE_PARALLEL_LOOP(MIN,MAX,CAPTURE(CAPTURES),INDEX,BODY)
  
#else
# define PAR_ON_EXEC_SPACE(EXEC_SPACE,MIN,MAX,CAPTURES,INDEX,BODY...)	\
  static_assert(EXEC_SPACE.hasUniqueExecSpace(),"Not uinique exec space!"); \
  HOST_PARALLEL_LOOP(MIN,MAX,CAPTURE(CAPTURES),INDEX,BODY)
#endif
}

#endif
