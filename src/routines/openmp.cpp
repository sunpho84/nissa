#include <omp.h>

#include "../base/debug.h"
#include "../base/global_variables.h"
#include "../base/init.h"
#include "../base/macros.h"
#include "../base/openmp_macros.h"
#include "ios.h"

 #include <unistd.h>

//put a barrier between threads
void thread_barrier(int barr_id)
{
#ifdef DEBUG
  if(IS_MASTER_THREAD) glb_barr_id=barr_id;
#endif

#pragma omp barrier

#ifdef DEBUG
  if(!IS_MASTER_THREAD)
    if(glb_barr_id!=barr_id) crash("Thread %d found barrier %d when waiting for %d",thread_id,barr_id,glb_barr_id);
#pragma omp barrier
#endif
}

//unlock the thread pool
void thread_pool_unlock()
{
  thread_barrier(UNLOCK_POOL_BARRIER);
#ifdef DEBUG
  printf("thread %d unlocking the pool\n",thread_id);
#endif
  thread_pool_locked=false;
}

//lock the thread pool
void thread_pool_lock()
{
  thread_pool_locked=true;
#ifdef DEBUG
  printf("thread %d locking the pool\n",thread_id);
#endif
  thread_barrier(LOCK_POOL_BARRIER);
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
void start_threaded_function(void(*function)(void))
{
#ifdef DEBUG
  printf("----------Start working thread pool--------\n");
#endif
  //set external function pointer and unlock pool threads
  threaded_function_ptr=function;
  thread_pool_unlock();
  
  //execute the function and relock the pool, so we are sure that they are not reading the work-to-do
  if(threaded_function_ptr!=NULL) threaded_function_ptr();
  thread_pool_lock();
  
#ifdef DEBUG
  printf("----------Finished working in the thread pool--------\n");
#endif
}

//delete the thread pool
void thread_pool_stop()
{
  //check to be thread 0
  if(thread_id!=0) crash("only thread 0 can stop the pool");
  
  //pass a NULL order
  start_threaded_function(NULL);
}

//start the master thread, locking all the other threads
void thread_master_start(int narg,char **arg,void(*main_function)(int narg,char **arg))
{
  //get the number of threads
  nthreads=omp_get_num_threads();
  master_printf("Starting %d threads\n",nthreads);
  
  //initialize nissa
  init_nissa(narg,arg);
  
  //launch the main function
  main_function(narg,arg);
  
  //exit the thread pool
  thread_pool_stop();
}

//start nissa in a threaded environment, sending all threads but first in the 
//thread pool and issuing the main function
void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg))
{
#pragma omp parallel
  {
    //take thread id and mark all thread but master as locked
    thread_id=omp_get_thread_num();
    thread_pool_locked=true;
    
    //distinguish master thread from the others
    if(thread_id!=0) thread_pool();
    else thread_master_start(narg,arg,main_function);
  }
}
