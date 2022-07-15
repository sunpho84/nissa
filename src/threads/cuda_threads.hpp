#ifndef _CUDA_THREADS_HPP
#define _CUDA_THREADS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <base/init.hpp>

#define NUM_THREADS 128

#define NACTIVE_THREADS 1
#define MANDATORY_PARALLEL
#define MANDATORY_NOT_PARALLEL
#define IS_PARALLEL false
#define GET_THREAD_ID()
#define THREAD_ID (0)
#define THREAD_BARRIER_FORCE()
#define THREAD_BARRIER()
#define IS_MASTER_THREAD (1)
// #define NISSA_PARALLEL_LOOP(INDEX,EXT_START,EXT_END) for(int64_t INDEX=EXT_START;INDEX<EXT_END;INDEX++)
// #define NISSA_PARALLEL_LOOP_END
#define NISSA_PARALLEL_LOOP_EXP(INDEX,EXT_START,EXT_END) cuda_parallel_for(__LINE__,__FILE__,EXT_START,EXT_END,[=] __device__ (const uint64_t& INDEX){
#define NISSA_PARALLEL_LOOP_END_EXP })
#define NISSA_PARALLEL_LOOP(INDEX,EXT_START,EXT_END) NISSA_PARALLEL_LOOP_EXP(INDEX,EXT_START,EXT_END)
#define NISSA_PARALLEL_LOOP_END NISSA_PARALLEL_LOOP_END_EXP
#define THREAD_ATOMIC_EXEC(inst) inst
#define THREAD_BROADCAST(out,in) (out)=(in)
#define THREAD_BROADCAST_PTR(out,in) THREAD_BROADCAST(out,in)

namespace nissa
{
  template <typename IMin,
	    typename IMax,
	    typename F>
  __global__
  void cuda_generic_kernel(const IMin min,
			   const IMax max,
			   F f)
  {
    const auto i=min+blockIdx.x*blockDim.x+threadIdx.x;
    if(i<max)
      f(i);
  }
  
  inline void thread_barrier_internal()
  {
    cudaDeviceSynchronize();
  }
  
  double take_time();
  
  template <typename IMin,
	    typename IMax,
	    typename F>
  void cuda_parallel_for(const int line,
			 const char *file,
			 const IMin min,
			 const IMax max,
			 F f)
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
	take_time();
      }
    
    if(length>0)
      {
	cuda_generic_kernel<<<grid_dimension,block_dimension>>>(min,max,std::forward<F>(f));
	thread_barrier_internal();
      }
    
    if(print)
      printf(" finished in %lg s\n",take_time()-initTime);
  }
  
  inline void cache_flush()
  {
  }
  
  inline double *glb_threads_reduce_double_vect(double *vect,int nel)
  {
    return vect;
  }
  
  //start nissa
  inline void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg),const char compile_info[5][1024])
  {
    //initialize nissa (master thread only)
    init_nissa(narg,arg,compile_info);
    
    main_function(narg,arg);
  }
}

#endif
