#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <omp.h>
#include <stdlib.h>

#define EXTERN_THREADS
 #include "openmp_threads.hpp"

#include "routines/ios.hpp"

#include "base/debug.hpp"
#include "base/init.hpp"
#include "base/vectors.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //start nissa in a threaded environment, sending all threads but first in the
  //thread pool and issuing the main function
  void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg),const char compile_info[5][1024])
  {
    //initialize nissa (master thread only)
    init_nissa(narg,arg,compile_info);
        
#pragma omp parallel
    {
      //get the number of threads and thread id
      nthreads=omp_get_num_threads();
      master_printf("Using %u threads\n",nthreads);
    }
    
    main_function(narg,arg);   
  }
}

