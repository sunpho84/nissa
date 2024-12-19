#ifndef _THREADS_HPP
#define _THREADS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <map>

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

  struct KernelBenchmark
  {
    std::string name;
    
    std::string file;
    
    int line;
    
    double totalTime;
    
    size_t nInvoke;
    
    KernelBenchmark(const std::string& name,
		    const std::string& file,
		    const int& line) :
      name(name),
      file(file),
      line(line),
      totalTime(0),
      nInvoke(0)
    {
    }
    
    KernelBenchmark() = default;
    
    KernelBenchmark(KernelBenchmark&&) = default;
  };
  
  inline std::map<size_t,KernelBenchmark> kernelBenchmarks;
  
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
    using I=
      std::common_type_t<IMin,IMax>;
    const I length=max-min;
    const dim3 blockDimension(128); // to be improved
    const dim3 gridDimension((length+blockDimension.x-1)/blockDimension.x);
    
    extern int rank;
    const bool print=(VERBOSITY_LV3
		      and rank==0);
    
    static const size_t hash=
      [line,
       file]()
      {
	const size_t hash=
	  hashCombine(typeid(F).hash_code(),line,std::string_view(file));
	
	const auto name=
	  typeid(F).name();
	
	kernelBenchmarks.try_emplace(hash,name,file,line);
	
	return hash;
      }();
    
    const double initTime=
      take_time();
    
    if(print)
      printf("at line %d of file %s launching kernel on loop [%d,%d) using blocks of size %d and grid of size %d\n",
	     line,file,(int)min,(int)max,blockDimension.x,gridDimension.x);
    
    if(length>0)
      {
	cudaGenericKernel<<<gridDimension,blockDimension>>>(min,max,f);
	cudaDeviceSynchronize();
      }
    
    const double runTime=
      take_time()-initTime;
    
    KernelBenchmark& b=
      kernelBenchmarks[hash];
    
    b.totalTime+=runTime;
    b.nInvoke++;
    
    if(print)
      printf(" finished in %lg s\n",runTime);
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
