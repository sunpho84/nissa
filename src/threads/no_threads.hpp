#ifndef _NO_THREADS_HPP
#define _NO_THREADS_HPP

#include "base/init.hpp"

#define NACTIVE_THREADS 1
#define MANDATORY_PARALLEL
#define MANDATORY_NOT_PARALLEL
#define IS_PARALLEL false
#define GET_THREAD_ID()
#define THREAD_ID (0)
#define THREAD_BARRIER_FORCE()
#define THREAD_BARRIER()
#define IS_MASTER_THREAD (1)
#define THREAD_ATOMIC_EXEC(inst) inst
#define THREAD_BROADCAST(out,in) (out)=(in)
#define THREAD_BROADCAST_PTR(out,in) THREAD_BROADCAST(out,in)

#define CUDA_MANAGED

#define CUDA_HOST_AND_DEVICE

#define NTHREADS 1

namespace nissa
{
  
  inline void cache_flush()
  {
  }
  
  inline double *glb_threads_reduce_double_vect(double *vect,int nel)
  {
    return vect;
  }
  
  inline void thread_barrier_internal()
  {
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
