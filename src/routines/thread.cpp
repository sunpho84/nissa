#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <omp.h>
#include <stdlib.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/init.hpp"
#include "base/macros.hpp"
#include "base/thread_macros.hpp"
#include "ios.hpp"

namespace nissa
{
  //put in the external bgq_barrier.c file, to avoid alignement problem
#if defined BGQ && (! defined BGQ_EMU)
#include "bgq/bgq_barrier.h"
#endif
  
  //put a barrier between threads
#ifdef THREAD_DEBUG
  void thread_barrier(const char *barr_file,int barr_line,int force_barrier)
#else
    void thread_barrier(int force_barrier)
#endif
  {
    if(!thread_pool_locked||force_barrier)
      {
	////////////////////////////////////// debug part 1 ///////////////////////////////////////////
#ifdef THREAD_DEBUG
	GET_THREAD_ID();
	if(nissa_verbosity>=3)
	  printf("thread %d rank %d barrier call on line %d of file %s (thread_pool_locked: %d, force_barrier: %d)\n",
		 thread_id,rank,barr_line,barr_file,thread_pool_locked,force_barrier);
	//copy the barrier id to the global ref
	if(IS_MASTER_THREAD)
	  {
	    strcpy(glb_barr_file,barr_file);
	    glb_barr_line=barr_line;
	    cache_flush();
	  }
#endif
	
	thread_barrier_internal();
	
	//////////////////////////////////// debug part 2 ///////////////////////////////////////////////
#ifdef THREAD_DEBUG
	//check that the local id correspond to global one
	if(!IS_MASTER_THREAD)
	  if(glb_barr_line!=barr_line||strcmp(glb_barr_file,barr_file))
	    crash("Thread %d found barrier on line %d of file %s when master thread invoked it at line %d of file %s)",
		  thread_id,barr_line,barr_file,glb_barr_line,glb_barr_file);
	thread_barrier_internal();
#endif
      }
  }
  
  //unlock the thread pool
  void thread_pool_unlock()
  {
    THREAD_BARRIER_FORCE();
#ifdef THREAD_DEBUG
    GET_THREAD_ID();
    if(rank==0) printf("thread %d unlocking the pool\n",thread_id);
#endif
    thread_pool_locked=false;
    cache_flush();
  }
  
  //lock the thread pool
  void thread_pool_lock()
  {
    THREAD_BARRIER_FORCE();
    thread_pool_locked=true;
    cache_flush();
#ifdef THREAD_DEBUG
    GET_THREAD_ID();
    if(rank==0) printf("thread %d locking the pool\n",thread_id);
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
    if(rank==0) printf("----------Start working %s thread pool--------\n",name);
#endif
    //set external function pointer and unlock pool threads
    threaded_function_ptr=function;
    thread_pool_unlock();
    
    //execute the function and relock the pool, so we are sure that they are not reading the work-to-do
    if(threaded_function_ptr!=NULL) threaded_function_ptr();
    thread_pool_lock();
    
#ifdef THREAD_DEBUG
    if(rank==0) printf("----------Finished working %s thread pool--------\n",name);
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
