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

#define MANDATORY_PARALLEL if(not IS_PARALLEL) crash("this cannot be called when threads are locked")
#define MANDATORY_NOT_PARALLEL if(IS_PARALLEL) crash("this cannot be called when threads are not locked")

#define IS_PARALLEL omp_in_parallel()

#define THREAD_ID omp_get_thread_num()
 
#define IS_MASTER_THREAD (THREAD_ID==0)

  template <typename T>
  auto __baseType(const T& t)
  {
    if constexpr(nissa::isTensorComp<T>)
      return t();
    else
      return t;
  }
  
template <typename A,typename B>
struct OMPThreadsIndices
{
  template <typename T>
  using BaseType=decltype(__baseType(std::declval<T>()));
  
  using ForType=std::common_type_t<BaseType<A>,BaseType<B>>;
  
  using ActualType=std::common_type_t<A,B>;
};

#define NISSA_PARALLEL_LOOP(INDEX,START,END)			\
  _Pragma("omp parallel for")					\
  for(typename OMPThreadsIndices<std::decay_t<decltype(START)>,std::decay_t<decltype(END)>>::ForType _ ## INDEX=__baseType(START);_ ## INDEX<__baseType(END);_ ## INDEX++){ \
  const typename OMPThreadsIndices<std::decay_t<decltype(START)>,std::decay_t<decltype(END)>>::ActualType INDEX= _ ## INDEX;
  
#define NISSA_PARALLEL_LOOP_END }
  
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
