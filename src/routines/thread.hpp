#ifndef _THREAD_HPP
#define _THREAD_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "new_types/float_128.hpp"
#include "base/random.hpp"

#ifndef EXTERN_THREAD
 #define EXTERN_THREAD extern
 #define INIT_TO(A)
 #define ONLY_INSTANTIATION
#else
#define INIT_TO(A) =A
#endif

namespace nissa
{
#ifdef THREAD_DEBUG
   EXTERN_THREAD int glb_barr_line;
   EXTERN_THREAD char glb_barr_file[1024];
  #if THREAD_DEBUG >=2
    EXTERN_THREAD rnd_gen *delay_rnd_gen;
    EXTERN_THREAD int *delayed_thread_barrier;
  #endif
 #endif
  
  EXTERN_THREAD bool thread_pool_locked INIT_TO(true);
  EXTERN_THREAD unsigned int nthreads INIT_TO(1);
  
  EXTERN_THREAD void *broadcast_ptr;
  EXTERN_THREAD float *glb_single_reduction_buf;
  EXTERN_THREAD double *glb_double_reduction_buf;
  EXTERN_THREAD float_128 *glb_quadruple_reduction_buf;
  
  EXTERN_THREAD void(*threaded_function_ptr)();
  
#ifdef THREAD_DEBUG
  void thread_barrier_with_check(const char*file,int line);
#else
  void thread_barrier_without_check();
#endif
  
  void start_threaded_function(void(*function)(void),const char *name);
  void thread_master_start(int narg,char **arg,void(*main_function)(int narg,char **arg));
  void thread_pool();
  void thread_pool_stop();
  double *glb_threads_reduce_double_vect(double *vect,int nel);
  inline complex *glb_threads_reduce_complex_vect(complex *vect,int nel)
  {return (complex*)glb_threads_reduce_double_vect((double*)vect,2*nel);}
}

#undef INIT_TO

#endif
