#ifndef _OPENMP_H
#define _OPENMP_H

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"

namespace nissa
{
#ifdef THREAD_DEBUG
  void thread_barrier_with_check(const char*file,int line);
#else
  void thread_barrier_without_check();
#endif
  
  void start_threaded_function(void(*function)(void),const char *name);
  void thread_master_start(int narg,char **arg,void(*main_function)(int narg,char **arg));
  void thread_pool();
  void thread_pool_stop();
  void glb_threads_reduce_double_vect(double *vect,int nel);
}

#endif
