#ifndef _THREADS_HPP
#define _THREADS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#ifdef USE_CUDA
# include <base/cuda.hpp>
#endif

#include <base/debug.hpp>

#include <metaprogramming/demangle.hpp>

#include <routines/mpi_routines.hpp>

#include <threads/benchmarks.hpp>

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
    using I=
      std::common_type_t<IMin,IMax>;
    
    double initTime=0;
    const bool print=
      (VERBOSITY_LV3
       and rank==0);
    
    if(print)
      {
	printf("at line %d of file %s launching openmp loop [%d,%d)\n",
	   line,file,(int)min,(int)max);
	initTime=take_time();
      }
    
#pragma omp parallel for
    for(I i=min;i<max;i++)
      f(IMax(i));
    
    if(print)
      printf(" finished in %lg s\n",take_time()-initTime);
  }
}

# define THREADS_PARALLEL_LOOP(ARGS...)		\
  openmpParallelFor(__LINE__,__FILE__,ARGS)

#else

# define THREADS_PARALLEL_LOOP(ARGS...) \
  SERIAL_LOOP(ARGS)

#endif

/////////////////////////////////////////////////////////////////

#ifdef USE_CUDA

namespace nissa
{
  /// Generic kernel to be used for the cuda launcher
  template <typename IMin,
	    typename IMax,
	    typename F>
  __global__
  void cudaGenericKernel(const IMin min,
			 const IMax max,
			 F f) // Needs to take by value
  {
    /// Type used for index, might we take it as int64_t?
    using I=
      std::common_type_t<IMin,IMax>;
    
    /// Index
    const I i=
      (min+blockIdx.x*blockDim.x+threadIdx.x);
    
    if(i<max)
      f(i);
  }
  
  /// Issue a parallel for using cuda kernel
  template <typename IMin,
	    typename IMax,
	    typename F>
  void cudaParallelFor(const int line,
		       const char* file,
		       const char* func,
		       const IMin min,
		       const IMax max,
		       F f) // Needs to take by value
  {
    /// Decide whether to print
    const bool print=
      (VERBOSITY_LV3
       and rank==0);
    
    /// Register the kernel and gets the id
    static const size_t id=
      [print,
       line,
       file,
       func]()
      {
	/// Name of the kernel
	const std::string name=
	  func+
	  (std::string)":"+
	  file+
	  (std::string)":"+
	  std::to_string(line);
	
	bool kernelIsFound=false;
	
	size_t id=0;
	while(id<kernelInfoLaunchParsStats.size() and not kernelIsFound)
	  if(kernelInfoLaunchParsStats[id].info.name==name)
	    kernelIsFound=true;
	  else
	    id++;
	
	if(not kernelIsFound)
	  {
	    kernelInfoLaunchParsStats.emplace_back(name);
	    if(print)
	      printf("Kernel %s not found, stored into id %zu/%zu\n",name.c_str(),id,kernelInfoLaunchParsStats.size());
	  }
	else
	  if(print)
	    printf("Found kernel %s, id: %zu\n",name.c_str(),id);
	
	return id;
      }();
    
    /// Type for the index of the loop
    using I=
      std::common_type_t<IMin,IMax>;
    
    /// Compute the length of the loop
    const I loopLength=
      max-min;
    
    if(loopLength>0)
      {
	VERBOSITY_LV3_MASTER_PRINTF("/////////////////////////////////////////////////////////////////\n");
	
	/// Attributes of the function
	static const cudaFuncAttributes attr=
	  []()
	  {
	    /// Fetch the attributes of the kernel function and return
	    cudaFuncAttributes attr;
	    cudaFuncGetAttributes(&attr,cudaGenericKernel<IMin,IMax,F>);
	    
	    return attr;
	  }();
	
	/// Proxy for the info of the kernel
	auto k=
	  []() -> KernelInfo&
	  {
	    return kernelInfoLaunchParsStats[id].info;
	  };
	
	VERBOSITY_LV3_MASTER_PRINTF("Going to launch kernel %s for loop size [%ld;%ld) = %ld\n",
		      k().name.c_str(),
		      (int64_t)min,
		      (int64_t)max,
		      (int64_t)loopLength);
	
	// printf("\n\n\nList of all kernel known and profiled\n");
	// for(const KernelInfoLaunchParsStat& kernelInfoLaunchParsStat : kernelInfoLaunchParsStats)
	//   {
	//     const auto& [info,launchParsStatList]=kernelInfoLaunchParsStat;
	//     for(auto& [loopLength,stat] : launchParsStatList)
	//       printf("name %s loopLength %ld optimal blocksize %d launched %ld times\n",
	// 	     info.name.c_str(),loopLength,stat.optimalBlockSize,stat.nInvoke);
	//   }
	// printf("\n\n\n\n");
	
	/// Launch the actual calculation
	auto launch=
	  [func,
	   file,
	   line,
	   min,
	   max,
	   loopLength,
	   f](const int blockSize)
	  {
	    /// Dimensions of the block
	    const dim3 blockDimension(blockSize);
	    
	    /// First entry in the sizes of the grid
	    const I g=(loopLength+blockDimension.x-1)/blockDimension.x;
	    
	    /// Size of the grid
	    const dim3 gridDimension(g);
	    
	    cudaGenericKernel<<<gridDimension,blockDimension>>>(min,max,f);
	    DECRYPT_CUDA_ERROR(cudaDeviceSynchronize(),"After launching kernel %s file %s line %d on loop [%d,%d) using blocks of size %d",
			       func,file,line,(int)min,(int)max,blockSize);
	  };
	
	/// Gets the maximum number of threads per block
	const int nMaxThreads=
	  attr.maxThreadsPerBlock;
	
	if(nMaxThreads<=0)
	  CRASH("kernel %s file %s line %d has max threads per block=%d",func,file,line,nMaxThreads);
	
	const int optimalBlockSize=
	  getOptimalBlockSize(id,loopLength,1,nMaxThreads,launch);
	
	if(optimalBlockSize<=0 or optimalBlockSize>nMaxThreads)
	  CRASH("Optimal number of threads %d is outside the allowed interval [0:%d]",optimalBlockSize,nMaxThreads);
	
	if(print)
	  printf("at line %d of file %s launching kernel on loop [%d,%d) using blocks of size %d\n",
		 line,file,(int)min,(int)max,optimalBlockSize);
	
	/// Take intiial time
	const double initTime=
	  take_time();
	
	VERBOSITY_LV3_MASTER_PRINTF("Launching with optimal blocksize %d\n",optimalBlockSize);
	launch(optimalBlockSize);
	
	/// Compute runtime
	const double runTime=
	  take_time()-initTime;
	
	KernelSizeLaunchParsStat& launchParsStat=
	  kernelInfoLaunchParsStats[id].launchParsStatList[loopLength];
	launchParsStat.totalTime+=runTime;
	launchParsStat.nInvoke++;
	
	VERBOSITY_LV3_MASTER_PRINTF(" finished in %lg s\n",runTime);
      }
  }
  
  /// Take notes of whether we are inside a parallel for
  ///
  /// Needed to take care of the reference
  inline int insideParallelFor;
}

# define CUDA_PARALLEL_LOOP(ARGS...)					\
  MACRO_GUARD(insideParallelFor++;					\
	      benchmarkBeginActions.clear();				\
	      benchmarkEndActions.clear();				\
	      cudaParallelFor(__LINE__,					\
			      __FILE__,					\
			      [](const auto& f)				\
			      {						\
				return typeid(f).name();		\
			      }([](){}),				\
			      ARGS);					\
	      insideParallelFor--;)

#else

# define CUDA_PARALLEL_LOOP(ARGS...)		\
  SERIAL_LOOP(ARGS)

#endif

#define _LAMBDA_FUNCTION_WITH_DECORATION_SUFFIX(DECORATION,		\
						SUFFIX,			\
						CAPTURES,		\
						INDEX,			\
						ARGS...)		\
  [CAPTURES] DECORATION (const auto& INDEX) SUFFIX			\
  ARGS
  
#define LAMBDA_FUNCTION_FOR_HOST(SUFFIX,				\
				 CAPTURES,				\
				 INDEX,					\
				 ARGS...)				\
  _LAMBDA_FUNCTION_WITH_DECORATION_SUFFIX(/* host */,			\
					  SUFFIX,			\
					  CAPTURE(CAPTURES),		\
					  INDEX,ARGS)

#define LAMBDA_FUNCTION_FOR_DEVICE(SUFFIX,				\
				   CAPTURES,				\
				   INDEX,				\
				   ARGS...)				\
  _LAMBDA_FUNCTION_WITH_DECORATION_SUFFIX(CUDA_DEVICE INLINE_ATTRIBUTE,	\
					  SUFFIX,			\
					  CAPTURE(CAPTURES),		\
					  INDEX,			\
					  ARGS)

#define HOST_PARALLEL_LOOP(MIN,						\
			   MAX,						\
			   CAPTURES,					\
			   INDEX,					\
			   ARGS...)					\
  THREADS_PARALLEL_LOOP(MIN,						\
			MAX,						\
			LAMBDA_FUNCTION_FOR_HOST(MUTABLE_INLINE_ATTRIBUTE, \
						 CAPTURE(CAPTURES),	\
						 INDEX,			\
						 ARGS))

#define DEVICE_PARALLEL_LOOP(MIN,					\
			     MAX,					\
			     CAPTURES,					\
			     INDEX,					\
			     ARGS...)					\
  CUDA_PARALLEL_LOOP(MIN,						\
		     MAX,						\
		     LAMBDA_FUNCTION_FOR_DEVICE(MUTABLE_INLINE_ATTRIBUTE, \
						CAPTURE(CAPTURES),	\
						INDEX,			\
						ARGS))

#if THREADS_TYPE == CUDA_THREADS

# define DEFAULT_PARALLEL_LOOP			\
  DEVICE_PARALLEL_LOOP

# define LAMBDA_FUNCTION_FOR_DEFAULT_THREADS	\
  LAMBDA_FUNCTION_FOR_DEVICE

#define LAMBDA_FUNCTION_DECORATION_FOR_DEFAULT_THREADS	\
  CUDA_DEVICE

#else

# define DEFAULT_PARALLEL_LOOP			\
  HOST_PARALLEL_LOOP

# define LAMBDA_FUNCTION_FOR_DEFAULT_THREADS	\
  LAMBDA_FUNCTION_FOR_HOST

#define LAMBDA_FUNCTION_DECORATION_FOR_DEFAULT_THREADS

#endif

#define PAR(MIN,					\
	    MAX,					\
	    CAPTURES,					\
	    ARGS...)					\
  DEFAULT_PARALLEL_LOOP(MIN,				\
			MAX,				\
			CAPTURE(CAPTURES),		\
			ARGS)

#endif
