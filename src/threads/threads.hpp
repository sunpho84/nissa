#ifndef _THREADS_HPP
#define _THREADS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#define NISSA_CHUNK_WORKLOAD(START,CHUNK_LOAD,END,EXT_START,EXT_END,CHUNK_ID,NCHUNKS) \
  int WORKLOAD=EXT_END-EXT_START,					\
  CHUNK_LOAD=(WORKLOAD+NCHUNKS-1)/NCHUNKS,				\
  START=EXT_START+CHUNK_ID*CHUNK_LOAD,					\
  END=START+CHUNK_LOAD< EXT_END ? START+CHUNK_LOAD : EXT_END

#define NISSA_CHUNK_LOOP(INDEX,EXT_START,EXT_END,CHUNK_ID,NCHUNKS)	\
  for(NISSA_CHUNK_WORKLOAD(START,CHUNK_LOAD,END,EXT_START,EXT_END,CHUNK_ID,NCHUNKS),INDEX=START;INDEX<END;INDEX++)


#if THREADS_TYPE == NO_THREADS
# include "threads/no_threads.hpp"
#elif THREADS_TYPE == OPENMP_THREADS
# include "threads/openmp_threads.hpp"
#elif THREADS_TYPE == CUDA_THREADS
# include "threads/cuda_threads.hpp"
#else
# error Unknown thread parallelization THREADS_TYPE !
#endif

#ifndef NISSA_PARALLEL_LOOP
 #define NISSA_PARALLEL_LOOP(INDEX,EXT_START,EXT_END) for(int INDEX=EXT_START;INDEX<EXT_END;INDEX++){
#endif

#ifndef NISSA_PARALLEL_LOOP_END
 #define NISSA_PARALLEL_LOOP_END }
#endif

#ifndef NISSA_PARALLEL_LOOP_EXP
 #define NISSA_PARALLEL_LOOP_EXP(INDEX,EXT_START,EXT_END) NISSA_PARALLEL_LOOP(INDEX,EXT_START,EXT_END)
#endif

#ifndef NISSA_PARALLEL_LOOP_END_EXP
 #define NISSA_PARALLEL_LOOP_END_EXP NISSA_PARALLEL_LOOP_END
#endif

/////////////////////////////////////////////////////////////////

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
  void openmp_parallel_for(const int line,
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

# define THREADS_PARALLEL_LOOP(ARGS...) openmp_parallel_for(__LINE__,__FILE__,ARGS)

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
  void cuda_generic_kernel(const IMin min,
			   const IMax max,
			   F& f)
  {
    const auto i=min+blockIdx.x*blockDim.x+threadIdx.x;
    if(i<max)
      f(i);
  }
  
  template <typename IMin,
	    typename IMax,
	    typename F>
  void cuda_parallel_for(const int line,
			 const char* file,
			 const IMin min,
			 const IMax max,
			 F&& f)
  {
    const auto length=(max-min);
    const dim3 block_dimension(NUM_THREADS);
    const dim3 grid_dimension((length+block_dimension.x-1)/block_dimension.x);
    
    double initTime=0;
    extern int rank,verbosity_lv;
    const bool print=(verbosity_lv>=1// 2
		      and rank==0);
    if(print)
      {
	printf("at line %d of file %s launching kernel on loop [%ld,%ld) using blocks of size %d and grid of size %d\n",
	   line,file,(int64_t)min,(int64_t)max,block_dimension.x,grid_dimension.x);
	initTime=take_time();
      }
    
    if(length>0)
      {
	cuda_generic_kernel<<<grid_dimension,block_dimension>>>(min,max,f);
	cudaDeviceSynchronize();
      }
    
    if(print)
      printf(" finished in %lg s\n",take_time()-initTime);
  }
}

# define CUDA_PARALLEL_LOOP(ARGS...) cuda_parallel_for(__LINE__,__FILE__,ARGS)

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
