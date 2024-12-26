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
    extern int rank;
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

# define THREADS_PARALLEL_LOOP(ARGS...) openmpParallelFor(__LINE__,__FILE__,ARGS)

#else

# define THREADS_PARALLEL_LOOP(ARGS...) SERIAL_LOOP(ARGS)

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
    using I=
      std::common_type_t<IMin,IMax>;
    
    const I i=
      (min+blockIdx.x*blockDim.x+threadIdx.x);
    
    if(i<max)
      f(i);
  }
  
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
	  func+(std::string)file+std::to_string(line);
	//typeid(F).name();
	
	bool kernelIdIsFound=false;
	
	size_t id=0;
	for(id=0;id<kernelInfoLaunchParsStats.size() and not kernelIdIsFound;id++)
	  if(kernelInfoLaunchParsStats[id].info.name==name)
	    kernelIdIsFound=true;
	
	if(not kernelIdIsFound)
	  {
	    kernelInfoLaunchParsStats.emplace_back(name,file,line);
	    if(print)
	      printf("Kernel %s not found, storing into id %zu, total size %zu\n",name.c_str(),id,kernelInfoLaunchParsStats.size());
	  }
	else
	  {
	    auto& i=
	      kernelInfoLaunchParsStats[id].info;
	    if(print)
	      printf("Providing kernel %s additional info on file %s and line %d\n",name.c_str(),file,line);
	    i.file=file;
	    i.line=line;
	  }
	
	return id;
      }();
    
    /// Attributes of the function
    static const cudaFuncAttributes attr=
      []()
      {
	/// Fetch the attributes of the kernel function and return
	cudaFuncAttributes attr;
	cudaFuncGetAttributes(&attr,cudaGenericKernel<IMin,IMax,F>);
	
	return attr;
      }();
    
    /// Type for the index of the loop
    using I=
      std::common_type_t<IMin,IMax>;
    
    /// Compute the length of the loop
    const I loopLength=
      max-min;
    
    /// Takes note of the action to carry out before benchmarking
    std::vector<BenchmarkAction> benchmarkBeginActions(std::move(nissa::benchmarkBeginActions));
    
    /// Takes note of the action to carry out after benchmarking
    std::vector<BenchmarkAction> benchmarkEndActions(std::move(nissa::benchmarkEndActions));
    
    if(loopLength>0)
      {
	master_printf("/////////////////////////////////////////////////////////////////\n");
	
	/// Proxy for the info of the kernel
	auto k=
	  []() -> KernelInfo&
	  {
	    return kernelInfoLaunchParsStats[id].info;
	  };
	
	master_printf("Going to launch kernel %s of file %s line %d for loop size [%ld;%ld) = %ld\n",
		      k().name.c_str(),
		      k().file.c_str(),
		      k().line,
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
	  [min,
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
	    cudaDeviceSynchronize();
	  };
	
	/// External value of the rank, needed to decide whether to print
	extern int rank;
	
	/// Gets the maximum number of threads per block
	const int nMaxThreads=
	  attr.maxThreadsPerBlock;
	
	const int optimalBlockSize=
	  getOptimalBlockSize(id,loopLength,1,nMaxThreads,launch);
	
	if(print)
	  printf("at line %d of file %s launching kernel on loop [%d,%d) using blocks of size %d\n",
		 line,file,(int)min,(int)max,optimalBlockSize);
	
	/// Take intiial time
	const double initTime=
	  take_time();
	
	master_printf("Launching with optimal blocksize %d\n",optimalBlockSize);
	launch(optimalBlockSize);
	
	/// Compute runtime
	const double runTime=
	  take_time()-initTime;
	
	KernelSizeLaunchParsStat& launchParsStat=
	  kernelInfoLaunchParsStats[id].launchParsStatList[loopLength];
	launchParsStat.totalTime+=runTime;
	launchParsStat.nInvoke++;
	
	if(print)
	  printf(" finished in %lg s\n",runTime);
      }
  }
  
  /// Take notes of whether we are inside a parallel for
  ///
  /// Needed to take care of the reference
  inline int insideParallelFor;
}

# define CUDA_PARALLEL_LOOP(ARGS...)					\
  MACRO_GUARD(insideParallelFor++;					\
	      cudaParallelFor(__LINE__,__FILE__,__FUNCTION__,ARGS);	\
	      insideParallelFor--;)

#else

# define CUDA_PARALLEL_LOOP(ARGS...) SERIAL_LOOP(ARGS)

#endif

#define HOST_PARALLEL_LOOP(MIN,MAX,CAPTURES,INDEX,A...) \
  THREADS_PARALLEL_LOOP(MIN,MAX,[CAPTURES] (const auto& INDEX) mutable A)

#define DEVICE_PARALLEL_LOOP(MIN,MAX,CAPTURES,INDEX,A...) \
  CUDA_PARALLEL_LOOP(MIN,MAX,[CAPTURES] CUDA_DEVICE INLINE_ATTRIBUTE (const auto& INDEX) mutable A)

#if THREADS_TYPE == CUDA_THREADS
# define DEFAULT_PARALLEL_LOOP DEVICE_PARALLEL_LOOP
#else
# define DEFAULT_PARALLEL_LOOP HOST_PARALLEL_LOOP
#endif
 
#define PAR(MIN,MAX,CAPTURES,A...)		\
  DEFAULT_PARALLEL_LOOP(MIN,MAX,CAPTURE(CAPTURES),A)

#endif
