#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <omp.h>
#include <stdlib.h>

#include "ios.hpp"

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/init.hpp"
#include "base/macros.hpp"
#include "base/random.hpp"
#include "base/thread_macros.hpp"
//put in the external bgq_barrier.c file, to avoid alignement problem
#if defined BGQ && (! defined BGQ_EMU)
 #include "bgq/bgq_barrier.hpp"
#endif

namespace nissa
{
#if THREAD_DEBUG>=2

  //wait previously delayed threads
  void wait_for_delayed_threads()
  {
    GET_THREAD_ID();
    
    if(delayed_thread_barrier[THREAD_ID]==0)
      {
	if(rank==0 && VERBOSITY_LV3) printf("thread %d waiting for delayed threads\n",THREAD_ID);
	thread_barrier_internal();
      }
  } 

  //select a new state for the delay
  void select_new_delay_pattern()
  {
    GET_THREAD_ID();
    
    //extract a random switch:
    // if 0, we exec immediately
    // if 1, we postpone to the end of the barrier
    delayed_thread_barrier[THREAD_ID]=(int)rnd_get_unif(delay_rnd_gen+THREAD_ID,0,2);
  }

  //delay marked threads
  void delay_marked_threads()
  {
    GET_THREAD_ID();
    
    if(delayed_thread_barrier[THREAD_ID]==1)
      {
	if(rank==0 && VERBOSITY_LV3) printf("thread %d will delay its execution\n",THREAD_ID);
	thread_barrier_internal();
      }
  }
#endif

#if THREAD_DEBUG>=1
  //check that the local barrier correspond to global one
  void check_barrier(const char *barr_file,int barr_line)
  {
    GET_THREAD_ID();
    
    if(VERBOSITY_LV3)
      printf("thread %d rank %d barrier call on line %d of file %s (thread_pool_locked: %d)\n",
	     thread_id,rank,barr_line,barr_file,thread_pool_locked);
    
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
	      thread_id,barr_line,barr_file,glb_barr_line,glb_barr_file);
  }
#endif

  //thread barrier without line number
#if THREAD_DEBUG>=1  
  void thread_barrier_with_check(const char *barr_file,int barr_line)
#else
  void thread_barrier_with_check()
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
  
  //unlock the thread pool
  void thread_pool_unlock()
  {
    THREAD_BARRIER_FORCE();
#ifdef THREAD_DEBUG
    GET_THREAD_ID();
    if(rank==0 && VERBOSITY_LV3) printf("thread %d unlocking the pool\n",thread_id);
#endif
    thread_pool_locked=false;
    cache_flush();
  }
  
  //lock the thread pool
  void thread_pool_lock()
  {
#ifdef THREAD_DEBUG
    GET_THREAD_ID();
#endif

    //if something was delayed, make it advance
    #if THREAD_DEBUG>=2
    wait_for_delayed_threads();
    delayed_thread_barrier[THREAD_ID]=0;
    #endif

    THREAD_BARRIER_FORCE();
    thread_pool_locked=true;
    cache_flush();
#ifdef THREAD_DEBUG
    if(rank==0 && VERBOSITY_LV3) printf("thread %d locking the pool\n",thread_id);
#endif
  }
  
  //thread pool (executed by non-master threads)
  void thread_pool()
  {
    GET_THREAD_ID();
    
    //check that thread 0 is not inside the pool
    if(thread_id==0) crash("thread 0 cannot enter the pool");
    
    //loop until asked to exit
    do
      {
	//hold until unlocked
	thread_pool_unlock();
	
	//wait that orders are communicated and in case exec them
	if(threaded_function_ptr!=NULL) threaded_function_ptr();

	thread_pool_lock();
      }
    while(threaded_function_ptr!=NULL);
    
#ifdef THREAD_DEBUG
    printf("thread %d exit pool\n",thread_id);
#endif
  }
  
  //execute a function using all threads
  void start_threaded_function(void(*function)(void),const char *name)
  {
#ifdef THREAD_DEBUG
    if(rank==0 && VERBOSITY_LV3) printf("----------Start working %s thread pool--------\n",name);
#endif
    //set external function pointer and unlock pool threads
    threaded_function_ptr=function;
    thread_pool_unlock();
    
    //execute the function and relock the pool, so we are sure that they are not reading the work-to-do
    if(threaded_function_ptr!=NULL) threaded_function_ptr();
    thread_pool_lock();
    
#ifdef THREAD_DEBUG
    if(rank==0 && VERBOSITY_LV3) printf("----------Finished working %s thread pool--------\n",name);
#endif
  }
  
  //delete the thread pool
  void thread_pool_stop()
  {
    GET_THREAD_ID();
    
    //check to be thread 0
    if(thread_id!=0) crash("only thread 0 can stop the pool");
    
    //pass a NULL order
    start_threaded_function(NULL,"");
  }
  
  //start the master thread, locking all the other threads
  void thread_master_start(int narg,char **arg,void(*main_function)(int narg,char **arg))
  {
    //initialize reducing buffers
    glb_double_reduction_buf=(double*)malloc(nthreads*sizeof(double));
    glb_float_128_reduction_buf=(float_128*)malloc(nthreads*sizeof(float_128));
    
    //launch the main function
    main_function(narg,arg);
    
    //free global reduction buffers
    free(glb_double_reduction_buf);
    free(glb_float_128_reduction_buf);

    //exit the thread pool
    thread_pool_stop();
  }
}
