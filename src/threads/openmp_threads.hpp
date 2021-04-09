#ifndef _OPENMP_THREADS_HPP
#define _OPENMP_THREADS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <omp.h>
#include <cstdint>

#include "base/debug.hpp"
#include "new_types/float_128.hpp"

#ifndef EXTERN_THREADS
 #define EXTERN_THREADS extern
 #define INIT_TO(A)
 #define ONLY_INSTANTIATION
#else
 #define INIT_TO(A) =A
#endif

#define CUDA_MANAGED

#define CUDA_HOST_DEVICE

#define NACTIVE_THREADS ((thread_pool_locked)?1:nthreads)
#define MANDATORY_PARALLEL if(nthreads>1 && thread_pool_locked) crash("this cannot be called when threads are locked")
#define MANDATORY_NOT_PARALLEL if(nthreads>1 && !thread_pool_locked) crash("this cannot be called when threads are not locked")

#define IS_PARALLEL (NACTIVE_THREADS!=1)

#define THREAD_ID omp_get_thread_num()
 
#ifdef THREAD_DEBUG
 #define THREAD_BARRIER_FORCE() thread_barrier_internal()
 #define THREAD_BARRIER()       if(!thread_pool_locked) thread_barrier_with_check(__FILE__,__LINE__)
#else
 #define THREAD_BARRIER_FORCE() thread_barrier_internal()
 #define THREAD_BARRIER()       if(!thread_pool_locked) thread_barrier_without_check()
#endif
 
#define IS_MASTER_THREAD (THREAD_ID==0)

#define NISSA_PARALLEL_LOOP(INDEX,START,END)			\
  NISSA_CHUNK_LOOP(INDEX,START,END,THREAD_ID,NACTIVE_THREADS){
#define NISSA_PARALLEL_LOOP_END }
  
#define THREAD_ATOMIC_EXEC(inst) do{THREAD_BARRIER();inst;THREAD_BARRIER();}while(0)
#define THREAD_BROADCAST(out,in)			\
  if(IS_MASTER_THREAD) broadcast_ptr=(void*)&in;		\
  THREAD_ATOMIC_EXEC(memcpy(&out,broadcast_ptr,sizeof(out)));
#define THREAD_BROADCAST_PTR(out,in)		\
  if(IS_MASTER_THREAD) broadcast_ptr=in;				\
  THREAD_ATOMIC_EXEC(memcpy(&out,&broadcast_ptr,sizeof(void*)));

namespace nissa
{
  //flush the cache
  inline void cache_flush()
  {
#pragma omp flush
  }
}

//barrier without any possible checking
namespace nissa
{
  inline void thread_barrier_internal()
  {
    #pragma omp barrier
  }
}

namespace nissa
{
#ifdef THREAD_DEBUG
   EXTERN_THREADS int glb_barr_line;
   EXTERN_THREADS char glb_barr_file[1024];
  #if THREAD_DEBUG >=2
    struct rnd_gen;
    EXTERN_THREADS rnd_gen *delay_rnd_gen;
    EXTERN_THREADS int *delayed_thread_barrier;
  #endif
 #endif
  
  EXTERN_THREADS bool thread_pool_locked INIT_TO(true);
  EXTERN_THREADS int nthreads INIT_TO(1);
  
  EXTERN_THREADS void *broadcast_ptr;
  
  EXTERN_THREADS void(*threaded_function_ptr)();
  
#ifdef THREAD_DEBUG
  void thread_barrier_with_check(const char*file,int line);
#else
  void thread_barrier_without_check();
#endif
  
  void start_threaded_function(void(*function)(void),const char *name);
  void thread_master_start(int narg,char **arg,void(*main_function)(int narg,char **arg));
  void thread_pool();
  void thread_pool_stop();
  
  void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg),const char compile_info[5][1024]);
}

#undef EXTERN_THREADS
#undef INIT_TO
#undef ONLY_INSTANTIATION

#endif
