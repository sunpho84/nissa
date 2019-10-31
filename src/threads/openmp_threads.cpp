#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <omp.h>
#include <stdlib.h>

#include "routines/ios.hpp"

#define EXTERN_THREADS
 #include "openmp_threads.hpp"

#include "base/debug.hpp"
#include "base/init.hpp"
#include "base/random.hpp"
#include "base/vectors.hpp"
#include "threads/threads.hpp"

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
    GET_THREAD_ID();
    
    if(delayed_thread_barrier[THREAD_ID]==1)
      {
	if(rank==0 && VERBOSITY_LV3) printf("thread %d will delay its execution,stopped at %s,%d\n",
					    THREAD_ID,glb_barr_file,glb_barr_line);
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
  
  //unlock the thread pool
  void thread_pool_unlock()
  {
    THREAD_BARRIER_FORCE();
#ifdef THREAD_DEBUG
    GET_THREAD_ID();
    if(rank==0 && VERBOSITY_LV3)
      {
	printf("thread %d unlocking the pool\n",THREAD_ID);
	fflush(stdout);
      }
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
    if(rank==0 && VERBOSITY_LV3)
      {
	printf("thread %d locking the pool\n",THREAD_ID);
	fflush(stdout);
      }
#endif
  }
  
  //thread pool (executed by non-master threads)
  void thread_pool()
  {
    GET_THREAD_ID();
    
    //check that thread 0 is not inside the pool
    if(THREAD_ID==0) crash("thread 0 cannot enter the pool");
    
    //set the thread as locked
    thread_pool_locked=true;
    
    //loop until asked to exit
    bool stay_working=true;
    do
      {
	//hold until unlocked
	thread_pool_unlock();
	
	//exec order or mark to exit in other case
	if(threaded_function_ptr!=NULL) threaded_function_ptr();
	else stay_working=false;
	
	thread_pool_lock();
      }
    while(stay_working);
    
#ifdef THREAD_DEBUG
    if(rank==0 && VERBOSITY_LV3)
      {
	printf("thread %d exit pool\n",THREAD_ID);
	fflush(stdout);
      }
#endif
  }
  
  //execute a function using all threads
  void start_threaded_function(void(*function)(void),const char *name)
  {
#ifdef THREAD_DEBUG
    if(rank==0 && VERBOSITY_LV3)
      {
	printf("----------Start working %s thread pool--------\n",name);
	fflush(stdout);
      }
#endif
    //set external function pointer and unlock pool threads
    threaded_function_ptr=function;
    thread_pool_unlock();
    
    //execute the function and relock the pool, so we are sure that they are not reading the work-to-do
    if(threaded_function_ptr!=NULL) threaded_function_ptr();
    thread_pool_lock();
    
#ifdef THREAD_DEBUG
    if(rank==0 && VERBOSITY_LV3)
      {
	printf("----------Finished working %s thread pool--------\n",name);
	fflush(stdout);
      }
#endif
  }
  
  //delete the thread pool
  void thread_pool_stop()
  {
    GET_THREAD_ID();
    
    //check to be thread 0
    if(THREAD_ID!=0) crash("only thread 0 can stop the pool");
    
    //pass a NULL order
    start_threaded_function(NULL,"");
  }
  
  //make all threads update a single counter in turn, checking previous state
  THREADABLE_FUNCTION_0ARG(thread_sanity_check)
  {
    //counter
    uint32_t *ptr=(uint32_t*)&broadcast_ptr;
    GET_THREAD_ID();
    
    //loop every threads, each one changing the counter state
    for(uint32_t i=0;i<NACTIVE_THREADS;i++)
      {
        //if it is current thread turn
        if(THREAD_ID==i)
          {
            //check that previous state agree (if not on first iter)
	    if(THREAD_ID!=0 && *ptr!=THREAD_ID-1)
              crash("error, thread %u found the counter in wrong state %u",THREAD_ID,*ptr);
            //update the counter
            *ptr=i;
          }
	
        THREAD_BARRIER();
      }
    
    //chech final state
    if(THREAD_ID==0 && *ptr!=nthreads-1) crash("loop thread not closed: %u!",*ptr);
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
  
  //start the master thread, locking all the other threads
  void thread_master_start(int narg,char **arg,void(*main_function)(int narg,char **arg))
  {
    //initialize reducing buffers
    glb_single_reduction_buf=(float*)malloc(nthreads*sizeof(float));
    glb_double_reduction_buf=(double*)malloc(nthreads*sizeof(double));
    glb_quadruple_reduction_buf=(float_128*)malloc(nthreads*sizeof(float_128));
    
    //lock the pool
    thread_pool_locked=true;
    cache_flush();
    
    //control the proper working of all the threads...
    thread_sanity_check();
    
    //launch the main function
    main_function(narg,arg);
    
    //free global reduction buffers
    free(glb_single_reduction_buf);
    free(glb_double_reduction_buf);
    free(glb_quadruple_reduction_buf);
    
    //exit the thread pool
    thread_pool_stop();
  }
  
  //reduce a double vector among threads
  double *glb_threads_reduce_double_vect(double *vect,int nel)
  {
    GET_THREAD_ID();
    
    //fill all the threads pointers
    double **ptr=nissa_malloc("ptr",NACTIVE_THREADS,double*);
    ptr[THREAD_ID]=vect;
    THREAD_BARRIER();
    
    //reduce among threads on thread 0
    for(unsigned int jthread=1;jthread<NACTIVE_THREADS;jthread++)
      {
	if(ptr[0]==ptr[jthread]) crash("overlapping vectors for threads %d and 0",jthread);
	NISSA_PARALLEL_LOOP(iel,0,nel)
	  ptr[0][iel]+=ptr[jthread][iel];
	NISSA_PARALLEL_LOOP_END;
      }
    
    //copy to the output
    for(unsigned int jthread=1;jthread<NACTIVE_THREADS;jthread++)
      NISSA_PARALLEL_LOOP(iel,0,nel)
	ptr[jthread][iel]=ptr[0][iel];
    NISSA_PARALLEL_LOOP_END;
    
    //return ptr 0
    double *ret=ptr[0];
    nissa_free(ptr);
    return ret;
  }
  
  //start nissa in a threaded environment, sending all threads but first in the
  //thread pool and issuing the main function
  void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg),const char compile_info[5][1024])
  {
    //initialize nissa (master thread only)
    init_nissa(narg,arg,compile_info);
    
    thread_pool_locked=false;
    cache_flush();
    
#pragma omp parallel
    {
      //get the number of threads and thread id
      nthreads=omp_get_num_threads();
      master_printf("Using %u threads\n",nthreads);
      
      //if BGQ, define appropriate barrier (to be done before sanity check, otherwise cannot barrier!)
#if defined BGQ && (! defined BGQ_EMU)
      bgq_barrier_define();
#endif
      
      //define delayed thread behavior (also this needed before sanity check, otherwise barrier would fail)
#if THREAD_DEBUG>=2
      delayed_thread_barrier=(int*)malloc(nthreads*sizeof(int));
      memset(delayed_thread_barrier,0,nthreads*sizeof(int));
      delay_rnd_gen=(rnd_gen*)malloc(nthreads*sizeof(rnd_gen));
      int delay_base_seed=time(0);
      for(unsigned int i=0;i<nthreads;i++) start_rnd_gen(delay_rnd_gen+i,delay_base_seed+i);
#endif
      
      //distinguish master thread from the others
      GET_THREAD_ID();
      if(thread_id!=0) thread_pool();
      else thread_master_start(narg,arg,main_function);
    }
  }
}
