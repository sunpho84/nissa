#ifndef _THREADS_HPP
#define _THREADS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <map>

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
  
  template <typename Head,
	    typename...Tail>
  inline size_t hashCombine(std::size_t seed,
			    const Head& v,
			    Tail&&...tail)
  {
    seed^=std::hash<Head>()(v)+0x9e3779b9+(seed<<6)+(seed>>2);
    
    if constexpr(sizeof...(Tail))
      return
	hashCombine(seed,std::forward<Tail>(tail)...);
    else
      return seed;
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
    /// Register the kernel and gets the id
    static const size_t id=
      [line,
       file]()
      {
	const std::string name=
	  demangle(typeid(F).name());
	
	kernelInfoLaunchParsStats.emplace_back(name,file,line);
	
	return kernelInfoLaunchParsStats.size()-1;
      }();
    
    for(const KernelInfoLaunchParsStat& kernelInfoLaunchParsStat : kernelInfoLaunchParsStats)
      {
	const auto& [info,launchParsStatList]=kernelInfoLaunchParsStat;
	for(auto& [size,stat] : launchParsStatList)
	  printf("name %s length %ld launched %ld times\n",info.name.c_str(),size,stat.nInvoke);
      }
    
    /// Type for the index
    using I=
      std::common_type_t<IMin,IMax>;
    
    /// Compute the length of the loop
    const I length=
      max-min;
    
    if(length>0)
      {
	/// Launch the actual calculation
	auto launch=
	  [&](const int blockSize)
	  {
	    const dim3 blockDimension(blockSize);
	    const dim3 gridDimension((length+blockSize-1)/blockSize);
	    
	    cudaGenericKernel<<<gridDimension,blockDimension>>>(min,max,f);
	    cudaDeviceSynchronize();
	  };
	
	extern int rank;
	
	const bool print=
	  (VERBOSITY_LV3
	   and rank==0);
	
	auto& [info,launchParsStatList]=
	  kernelInfoLaunchParsStats[id];
	
	auto [launchParsStatIt,toBeComputed]=
	  launchParsStatList.try_emplace(length);
	
	/// Reference to the parameters for the launch
	KernelSizeLaunchParsStat& launchParsStat=
	  launchParsStatIt->second;
	
	/// Reference to the optimal block size
	int& optimalBlockSize=
	  launchParsStat.optimalBlockSize;
	
	printf("ToBeComputed: %d\n",toBeComputed);
	
	if(toBeComputed)
	  {
	    printf("Benchmarking kernel %s for length %ld\n",
		   info.name.c_str(),(int64_t)length);
	    
	    toBeComputed=false;
	    
	    benchmarkInProgress=true;
	    
	    for(BenchmarkAction& b : benchmarkBeginActions)
	      b();
	    benchmarkBeginActions.clear();
	    
	    optimalBlockSize=0;
	    
	    const int nBench=100;
	    double minTime=0.0;
	    
	    printf("starting testBlockSize\n");
	    for(int testBlockSize=64;testBlockSize<=1024;testBlockSize*=2)
	      {
		// warmup
		launch(testBlockSize);
		
		const double initTime=
		  take_time();
		
		for(int i=0;i<nBench;i++)
		  launch(testBlockSize);
		
		const double runTime=
		  take_time()-initTime;
		
		if(optimalBlockSize==0 or minTime>runTime)
		  {
		    optimalBlockSize=testBlockSize;
		    minTime=runTime;
		  }
		
		printf("Benchmarked with blockSize %d, runtime %lg s minimal %lg s current optimal size %d\n",testBlockSize,runTime,minTime,optimalBlockSize);
	      }
	    
	    benchmarkInProgress=false;
	    
	    for(BenchmarkAction& e : benchmarkEndActions)
	      e();
	    benchmarkEndActions.clear();
	  }
	
	if(print)
	  printf("at line %d of file %s launching kernel on loop [%d,%d) using blocks of size %d\n",
		 line,file,(int)min,(int)max,optimalBlockSize);
	
	/// Take intiial time
	const double initTime=
	  take_time();
	
	launch(optimalBlockSize);
	
	/// Compute runtime
	const double runTime=
	  take_time()-initTime;
	
	launchParsStat.totalTime+=runTime;
	launchParsStat.nInvoke++;
	
	if(print)
	  printf(" finished in %lg s\n",runTime);
      }
  }
}

# define CUDA_PARALLEL_LOOP(ARGS...) cudaParallelFor(__LINE__,__FILE__,ARGS)

#else

# define CUDA_PARALLEL_LOOP(ARGS...) SERIAL_LOOP(ARGS)

#endif

#define HOST_PARALLEL_LOOP(MIN,MAX,CAPTURES,INDEX,A...) \
  THREADS_PARALLEL_LOOP(MIN,MAX,[CAPTURES] (const auto& INDEX) mutable A)

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
