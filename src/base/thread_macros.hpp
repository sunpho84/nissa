#ifndef _THREAD_MACROS_HPP
#define _THREAD_MACROS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "debug.hpp"
#ifdef USE_THREADS
 #include <omp.h>
 #include "routines/thread.hpp"
#endif

//////////////////////////////////////////////////////////////////////////////////////

#ifndef USE_THREADS
 #define NACTIVE_THREADS 1
 #define MANDATORY_PARALLEL
 #define MANDATORY_NOT_PARALLEL
#else
 #define NACTIVE_THREADS ((thread_pool_locked)?1:nthreads)
 #define MANDATORY_PARALLEL if(nthreads>1 && thread_pool_locked) CRASH("this cannot be called when threads are locked")
 #define MANDATORY_NOT_PARALLEL if(nthreads>1 && !thread_pool_locked) CRASH("this cannot be called when threads are not locked")
#endif

//barriers and flush
#ifdef USE_THREADS
 #define cache_flush() _Pragma("omp flush")
 #define thread_barrier_internal() _Pragma("omp barrier")
#else
 #define cache_flush()
 #define thread_barrier_internal()
#endif

#define IS_PARALLEL (NACTIVE_THREADS!=1)

#define NISSA_CHUNK_WORKLOAD(START,CHUNK_LOAD,END,EXT_START,EXT_END,CHUNK_ID,NCHUNKS) \
  int WORKLOAD=EXT_END-EXT_START,					\
  CHUNK_LOAD=(WORKLOAD+NCHUNKS-1)/NCHUNKS,				\
  START=EXT_START+CHUNK_ID*CHUNK_LOAD,					\
  END=START+CHUNK_LOAD< EXT_END ? START+CHUNK_LOAD : EXT_END

#define NISSA_CHUNK_LOOP(INDEX,EXT_START,EXT_END,CHUNK_ID,NCHUNKS)	\
  for(NISSA_CHUNK_WORKLOAD(START,CHUNK_LOAD,END,EXT_START,EXT_END,CHUNK_ID,NCHUNKS),INDEX=START;INDEX<END;INDEX++)

#ifdef USE_THREADS

 #define GET_THREAD_ID() uint32_t thread_id=omp_get_thread_num()
 #define THREAD_ID thread_id
 
 #ifdef THREAD_DEBUG
  #define THREAD_BARRIER_FORCE() thread_barrier_internal()
  #define THREAD_BARRIER()       if(!thread_pool_locked) thread_barrier_with_check(__FILE__,__LINE__)
 #else
  #define THREAD_BARRIER_FORCE() thread_barrier_internal()
  #define THREAD_BARRIER()       if(!thread_pool_locked) thread_barrier_without_check()
 #endif
 
 #define IS_MASTER_THREAD (THREAD_ID==0)
 
 #define NISSA_PARALLEL_LOOP(INDEX,START,END)			\
  NISSA_CHUNK_LOOP(INDEX,START,END,thread_id,NACTIVE_THREADS)
 
 #define THREAD_ATOMIC_EXEC(inst) do{THREAD_BARRIER();inst;THREAD_BARRIER();}while(0)
 #define THREAD_BROADCAST(out,in)			\
   if(IS_MASTER_THREAD) broadcast_ptr=(void*)&in;	\
   THREAD_ATOMIC_EXEC(memcpy(&out,broadcast_ptr,sizeof(out)));
 #define THREAD_BROADCAST_PTR(out,in)		\
   if(IS_MASTER_THREAD) broadcast_ptr=in;	\
   THREAD_ATOMIC_EXEC(memcpy(&out,&broadcast_ptr,sizeof(void*)));
  
#else

 #define GET_THREAD_ID()
 #define THREAD_ID (0)
 #define THREAD_BARRIER_FORCE()
 #define THREAD_BARRIER()
 #define IS_MASTER_THREAD (1)
 #define NISSA_PARALLEL_LOOP(INDEX,EXT_START,EXT_END) for(int INDEX=EXT_START;INDEX<EXT_END;INDEX++)
 #define THREAD_ATOMIC_EXEC(inst) inst
 #define THREAD_BROADCAST(out,in) (out)=(in)
 #define THREAD_BROADCAST_PTR(out,in) THREAD_BROADCAST(out,in)

#endif

////////////////////// HEADERS FOR THREADED FUNCTIONS ///////////////////////////////

#define THREADABLE_FUNCTION_0ARG_HEADER(FUNC_NAME) void FUNC_NAME()
#define THREADABLE_FUNCTION_1ARG_HEADER(FUNC_NAME,AT1,A1) void FUNC_NAME(AT1 A1)
#define THREADABLE_FUNCTION_2ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2) void FUNC_NAME(AT1 A1,AT2 A2)
#define THREADABLE_FUNCTION_3ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3)
#define THREADABLE_FUNCTION_4ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4)
#define THREADABLE_FUNCTION_5ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5)
#define THREADABLE_FUNCTION_6ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6)
#define THREADABLE_FUNCTION_7ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7)
#define THREADABLE_FUNCTION_8ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8)
#define THREADABLE_FUNCTION_9ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8,AT9 A9)
#define THREADABLE_FUNCTION_10ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8,AT9 A9,AT10 A10)
#define THREADABLE_FUNCTION_11ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8,AT9 A9,AT10 A10,AT11 A11)

//////////////////////////////////////////////////////////////////////////////////////

#ifdef USE_THREADS

//thread debug
#ifdef THREAD_DEBUG
 #define DEBUG_THREAD_START(NAME)					\
  do									\
    if(rank==0 && VERBOSITY_LV3)					\
      {									\
	printf("----------Start working %s thread pool--------\n",NAME); \
	fflush(stdout);							\
      }									\
  while(0)
 #define DEBUG_THREAD_END(NAME)						\
  do									\
    if(rank==0 && VERBOSITY_LV3)					\
      {									\
	printf("----------Finished working %s thread pool--------\n",NAME); \
	fflush(stdout);							\
      }									\
  while(0)
#else
 #define DEBUG_THREAD_START(NAME)
 #define DEBUG_THREAD_END(NAME)
#endif

#define THREADABLE_FUNCTION_SUMMON(FUNC_NAME)				\
  {									\
    if(nthreads>1 && thread_pool_locked)				\
      {									\
	DEBUG_THREAD_START(#FUNC_NAME);					\
	threaded_function=[&](){

#define ED(FUNC_NAME)							\
  ;};									\
  thread_pool_unlock();						        \
  threaded_function();							\
  thread_pool_lock();							\
  DEBUG_THREAD_END(FUNC_NAME)						\
  }									\
  else

#else

//skip execution
#define THREADABLE_FUNCTION_SUMMON(FUNC_NAME) { if(0)
#define ED ;

#endif

#define THREADABLE_FUNCTION_0ARG(FUNC_NAME) \
  THREADABLE_FUNCTION_0ARG_HEADER(FUNC_NAME)	    \
  THREADABLE_FUNCTION_SUMMON(FUNC_NAME) FUNC_NAME() ED(FUNC_NAME)
#define THREADABLE_FUNCTION_1ARG(FUNC_NAME,AT1,A1) \
  THREADABLE_FUNCTION_1ARG_HEADER(FUNC_NAME,AT1,A1) \
  THREADABLE_FUNCTION_SUMMON(FUNC_NAME) FUNC_NAME(A1) ED(FUNC_NAME)
#define THREADABLE_FUNCTION_2ARG(FUNC_NAME,AT1,A1,AT2,A2) \
  THREADABLE_FUNCTION_2ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2) \
  THREADABLE_FUNCTION_SUMMON(FUNC_NAME) FUNC_NAME(A1,A2) ED(FUNC_NAME)
#define THREADABLE_FUNCTION_3ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3) \
  THREADABLE_FUNCTION_3ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3) \
  THREADABLE_FUNCTION_SUMMON(FUNC_NAME) FUNC_NAME(A1,A2,A3) ED(FUNC_NAME)
#define THREADABLE_FUNCTION_4ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4) \
  THREADABLE_FUNCTION_4ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4) \
  THREADABLE_FUNCTION_SUMMON(FUNC_NAME) FUNC_NAME(A1,A2,A3,A4) ED(FUNC_NAME)
#define THREADABLE_FUNCTION_5ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) \
  THREADABLE_FUNCTION_5ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) \
  THREADABLE_FUNCTION_SUMMON(FUNC_NAME) FUNC_NAME(A1,A2,A3,A4,A5) ED(FUNC_NAME)
#define THREADABLE_FUNCTION_6ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) \
  THREADABLE_FUNCTION_6ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) \
  THREADABLE_FUNCTION_SUMMON(FUNC_NAME) FUNC_NAME(A1,A2,A3,A4,A5,A6) ED(FUNC_NAME)
#define THREADABLE_FUNCTION_7ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) \
  THREADABLE_FUNCTION_7ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) \
  THREADABLE_FUNCTION_SUMMON(FUNC_NAME) FUNC_NAME(A1,A2,A3,A4,A5,A6,A7) ED(FUNC_NAME)
#define THREADABLE_FUNCTION_8ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) \
  THREADABLE_FUNCTION_8ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) \
  THREADABLE_FUNCTION_SUMMON(FUNC_NAME) FUNC_NAME(A1,A2,A3,A4,A5,A6,A7,A8) ED(FUNC_NAME)
#define THREADABLE_FUNCTION_9ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) \
  THREADABLE_FUNCTION_9ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) \
  THREADABLE_FUNCTION_SUMMON(FUNC_NAME) FUNC_NAME(A1,A2,A3,A4,A5,A6,A7,A8,A9) ED(FUNC_NAME)
#define THREADABLE_FUNCTION_10ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) \
  THREADABLE_FUNCTION_10ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) \
  THREADABLE_FUNCTION_SUMMON(FUNC_NAME) FUNC_NAME(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10) ED(FUNC_NAME)
#define THREADABLE_FUNCTION_11ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11) \
  THREADABLE_FUNCTION_11ARG_HEADER(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11) \
  THREADABLE_FUNCTION_SUMMON(FUNC_NAME) FUNC_NAME(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11) ED(FUNC_NAME)

#define THREADABLE_FUNCTION_END }

#define FORM_TWO_THREAD_TEAMS()						\
  bool is_in_first_team,is_in_second_team;				\
  unsigned int nthreads_in_team,thread_in_team_id;			\
  if(thread_pool_locked||nthreads==1)					\
    {									\
      is_in_first_team=is_in_second_team=true;				\
      nthreads_in_team=1;						\
      thread_in_team_id=0;						\
    }									\
  else									\
    {									\
      is_in_first_team=(thread_id<nthreads/2);				\
      is_in_second_team=!is_in_first_team;				\
      if(is_in_first_team)						\
	{								\
	  nthreads_in_team=nthreads/2;					\
	  thread_in_team_id=thread_id;					\
	}								\
      else								\
	{								\
	  nthreads_in_team=nthreads-nthreads/2;				\
	  thread_in_team_id=thread_id-nthreads/2;			\
	}								\
    }
#endif
