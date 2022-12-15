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

#define NTHREADS nthreads

EXTERN_THREADS int nthreads INIT_TO(1);

#define CUDA_MANAGED

#define CUDA_HOST_AND_DEVICE

#define MANDATORY_PARALLEL if(not IS_PARALLEL) crash("this cannot be called when threads are locked")
#define MANDATORY_NOT_PARALLEL if(IS_PARALLEL) crash("this cannot be called when threads are not locked")

#define IS_PARALLEL omp_in_parallel()

#define THREAD_ID omp_get_thread_num()
 
#define IS_MASTER_THREAD (THREAD_ID==0)

#define NISSA_PARALLEL_LOOP(INDEX,START,END)			\
  _Pragma("omp parallel for")					\
  for(std::common_type_t<std::decay_t<decltype(START)>,std::decay_t<decltype(END)>> INDEX=START;INDEX<END;INDEX++){
#define NISSA_PARALLEL_LOOP_END }

namespace nissa
{
  template <typename IMin,
	    typename IMax,
	    typename F>
  void openmp_parallel_for(const int line,
			   const char *file,
			   const IMin min,
			   const IMax max,
			   F f)
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

#define PARALLEL_LOOP(ARGS...) openmp_parallel_for(__LINE__,__FILE__,ARGS)

#define THREAD_ATOMIC_EXEC(inst) do{THREAD_BARRIER();inst;THREAD_BARRIER();}while(0)

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

# define THREAD_BARRIER()      \
  _Pragma("omp barrier")

namespace nissa
{
  void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg),const char compile_info[5][1024]);
}

#undef EXTERN_THREADS
#undef INIT_TO
#undef ONLY_INSTANTIATION

#endif
