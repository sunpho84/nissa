#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <omp.h>
#include <stdlib.h>

#include "../base/debug.h"
#include "../base/global_variables.h"
#include "../base/init.h"
#include "../base/macros.h"
#include "../base/openmp_macros.h"
#include "ios.h"

//#define DEBUG

//put in the external bgq_barrier.c file, to avoid alignement problem
#ifdef BGQ
extern "C" void bgq_barrier(int n);
extern "C" void bgq_barrier_define();
#endif

//put a barrier between threads
void thread_barrier(int barr_id,int force_barrier=false)
{

  if(!thread_pool_locked||force_barrier)
    {
#ifdef DEBUG
      //debug: copy the barrier id to the global ref
      if(IS_MASTER_THREAD)
	{
	  glb_barr_id=barr_id;
#pragma omp flush(glb_barr_id)
	}
#endif
      
      //barrier
#ifdef BGQ
      bgq_barrier(nthreads);
#else
      #pragma omp barrier
#endif
      
#ifdef DEBUG
      //debug: check that the local id correspond to global one
      if(!IS_MASTER_THREAD)
	if(glb_barr_id!=barr_id) crash("Thread %d found barrier %d when waiting for %d",thread_id,barr_id,glb_barr_id);
#ifdef BGQ
      bgq_barrier(nthreads);
#else
      #pragma omp barrier
#endif
#endif
    }
}

//unlock the thread pool
void thread_pool_unlock()
{
  THREAD_BARRIER_FORCE(UNLOCK_POOL_BARRIER);
#ifdef DEBUG
  if(rank==0) printf("thread %d unlocking the pool\n",thread_id);
#endif
  thread_pool_locked=false;
#pragma omp flush(thread_pool_locked)
}

//lock the thread pool
void thread_pool_lock()
{
  THREAD_BARRIER_FORCE(LOCK_POOL_BARRIER);
  thread_pool_locked=true;
#pragma omp flush(thread_pool_locked)
#ifdef DEBUG
  if(rank==0) printf("thread %d locking the pool\n",thread_id);
#endif
}

//thread pool (executed by non-master threads)
void thread_pool()
{
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
  
#ifdef DEBUG
  printf("thread %d exit pool\n",thread_id);
#endif
}

//execute a function using all threads
void start_threaded_function(void(*function)(void),char *name)
{
#ifdef DEBUG
  if(rank==0) printf("----------Start working %s thread pool--------\n",name);
#endif
  //set external function pointer and unlock pool threads
  threaded_function_ptr=function;
  thread_pool_unlock();
  
  //execute the function and relock the pool, so we are sure that they are not reading the work-to-do
  if(threaded_function_ptr!=NULL) threaded_function_ptr();
  thread_pool_lock();
  
#ifdef DEBUG
  if(rank==0) printf("----------Finished working in the thread pool--------\n");
#endif
}

//delete the thread pool
void thread_pool_stop()
{
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

//start nissa in a threaded environment, sending all threads but first in the 
//thread pool and issuing the main function
void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg))
{
  //if BGQ, define appropriate barrier
#ifdef BGQ
  bgq_barrier_define();
#endif

#pragma omp parallel
  {
    //initialize nissa (master thread only)
#pragma omp master
    init_nissa(narg,arg);
#pragma omp barrier
  
    //get the number of threads and thread id
    nthreads=omp_get_num_threads();
    master_printf("Using %d threads\n",nthreads);
    
    //distinguish master thread from the others
    if(thread_id!=0) thread_pool();
    else thread_master_start(narg,arg,main_function);
  }
}
