#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <omp.h>
#include <stdlib.h>

#define EXTERN_THREADS
 #include "openmp_threads.hpp"

#include "routines/ios.hpp"

#include "base/debug.hpp"
#include "base/init.hpp"
#include "base/random.hpp"
#include "base/vectors.hpp"
#include "threads/threads.hpp"

namespace nissa
{
#if THREAD_DEBUG>=2
  
  //wait previously delayed threads
  void wait_for_delayed_threads()
  {
    
    if(delayed_thread_barrier[THREAD_ID]==0)
      {
	if(is_master_rank() && VERBOSITY_LV3) printf("thread %d waiting for delayed threads\n",THREAD_ID);
	thread_barrier_internal();
      }
  }
  
  //select a new state for the delay
  void select_new_delay_pattern()
  {
    
    //extract a random switch:
    // if 0, we exec immediately
    // if 1, we postpone to the end of the barrier
    
    enum delay_pattern{DELAY_RANDOMLY,DELAY_SLAVES};
    const delay_pattern picked=DELAY_RANDOMLY;
    
    switch(picked)
      {
      case DELAY_RANDOMLY:
	delayed_thread_barrier[THREAD_ID]=(int)rnd_get_unif(delay_rnd_gen+THREAD_ID,0,2);
	break;
      case DELAY_SLAVES:
	delayed_thread_barrier[THREAD_ID]=!IS_MASTER_THREAD;
	break;
      default:
	crash("Unknown delay pattern %d",picked);
      }
  }
  
  //delay marked threads
  void delay_marked_threads()
  {
    
    if(delayed_thread_barrier[THREAD_ID]==1)
      {
	if(is_master_rank() && VERBOSITY_LV3) printf("thread %d will delay its execution,stopped at %s,%d\n",
					    THREAD_ID,glb_barr_file,glb_barr_line);
	thread_barrier_internal();
      }
  }
#endif
  
#if THREAD_DEBUG>=1
  //check that the local barrier correspond to global one
  void check_barrier(const char *barr_file,int barr_line)
  {
    
    if(VERBOSITY_LV3)
      printf("thread %d rank %d barrier call on line %d of file %s (thread_pool_locked: %d)\n",
	     THREAD_ID,rank,barr_line,barr_file,thread_pool_locked);
    
    //copy the barrier id to the global ref
    if(IS_MASTER_THREAD)
      {
	strcpy(glb_barr_file,barr_file);
	glb_barr_line=barr_line;
	cache_flush();
      }
    thread_barrier_internal();
    
    //check
    if(!IS_MASTER_THREAD)
      if(glb_barr_line!=barr_line||strcmp(glb_barr_file,barr_file))
	crash("Thread %d found barrier on line %d of file %s when master thread invoked it at line %d of file %s)",
	      THREAD_ID,barr_line,barr_file,glb_barr_line,glb_barr_file);
  }
#endif

  //thread barrier without line number
#if THREAD_DEBUG>=1
  void thread_barrier_with_check(const char *barr_file,int barr_line)
#else
  void thread_barrier_without_check()
#endif
  {
    //if something was delayed, make it advance
    #if THREAD_DEBUG>=2
    wait_for_delayed_threads();
    #endif
    
    //true barrier
    thread_barrier_internal();
    
    //check coherence of called barrier
    #if THREAD_DEBUG>=1
    check_barrier(barr_file,barr_line);
    #endif
    
    //selct new delay pattern and delay
    #if THREAD_DEBUG>=2
    select_new_delay_pattern();
    delay_marked_threads();
    #endif
  }
  
  //make all threads update a single counter in turn, checking previous state
  //start nissa in a threaded environment, sending all threads but first in the
  //thread pool and issuing the main function
  void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg),const char compile_info[5][1024])
  {
    //initialize nissa (master thread only)
    init_nissa(narg,arg,compile_info);
    
    cache_flush();
    
#pragma omp parallel
    {
      //get the number of threads and thread id
      nthreads=omp_get_num_threads();
      master_printf("Using %u threads\n",nthreads);
      
      //define delayed thread behavior (also this needed before sanity check, otherwise barrier would fail)
#if THREAD_DEBUG>=2
      delayed_thread_barrier=(int*)malloc(nthreads*sizeof(int));
      memset(delayed_thread_barrier,0,nthreads*sizeof(int));
      delay_rnd_gen=(rnd_gen*)malloc(nthreads*sizeof(rnd_gen));
      int delay_base_seed=time(0);
      for(unsigned int i=0;i<nthreads;i++) start_rnd_gen(delay_rnd_gen+i,delay_base_seed+i);
#endif
      
      //distinguish master thread from the others
      main_function(narg,arg);
    }
  }
}

