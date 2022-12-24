#ifndef _CUDA_THREADS_HPP
#define _CUDA_THREADS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <base/init.hpp>

#define NUM_THREADS 128

#define NTHREADS 1
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

namespace nissa
{
  inline void thread_barrier_internal()
  {
    cudaDeviceSynchronize();
  }
  
  double take_time();
  
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
