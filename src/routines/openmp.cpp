#include <omp.h>

#include "../base/debug.h"
#include "../base/global_variables.h"
#include "../base/init.h"
#include "../base/macros.h"
#include "ios.h"

//put a barrier between threads
void thread_barrier()
{
#pragma omp barrier
}

//thread pool
void thread_pool()
{
  //check that thread 0 is not inside the ppol
  if(thread_id==0) crash("thread 0 must not enter the pool");
  
  //loop until not asked to exit
  do
    {
      //hold until not released, then if we have to stay in pool execute function
      thread_barrier();
      if(thread_stay_in_pool)
	{
	  threaded_function_ptr();
	}
    }
  while(thread_stay_in_pool);
}

//execute a function using all threads
void threaded_function_exec(void(*function)(void))
{
  //set external function pointer and release threads
  threaded_function_ptr=function;
  thread_barrier();
  
  //execute the function
  threaded_function_ptr();
}

//delete the thread pool
void thread_pool_stop()
{
  //check to be thread 0
  if(thread_id!=0) crash("only thread 0 can stop the pool");
  
  //release threads from the pool
  thread_stay_in_pool=0;
  
  //now we just have to release the other threads
  thread_barrier();
}

//start the master thread, keeping all the other threads in hold
void thread_master_start(int narg,char **arg,void(*main_function)(int narg,char **arg))
{
  //write the number of thread
  master_printf("Starting %d threads\n",omp_get_num_threads());
  
  //stuck the other thread in the pool
  thread_stay_in_pool=1;
  
  //initialize nissa
  init_nissa(narg,arg);
  
  //launch the main function
  main_function(narg,arg);
  
  //exit from the thread pool
  thread_pool_stop();
}

//start nissa in a threaded environment, sending all threads but first in the 
//thread pool and issuing the main function
void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg))
{
#pragma omp parallel
  {
    //take local thread id
    thread_id=omp_get_thread_num();
    
    //discriminate master thread from the others
    if(thread_id!=0) thread_pool();
    else thread_master_start(narg,arg,main_function);
  }
}
