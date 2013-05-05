#ifndef _OPENMP_MACROS_H
#define _OPENMP_MACROS_H

#include "debug.h"
#include <omp.h>

///////////////////////////////////// barriers ////////////////////////////////////////

#define LOCK_POOL_BARRIER 1
#define UNLOCK_POOL_BARRIER 2
#define SET_VEC_FLAG_FIRST_BARRIER 3
#define SET_VEC_FLAG_SECOND_BARRIER 4
#define UNSET_VEC_FLAG_FIRST_BARRIER 5
#define UNSET_VEC_FLAG_SECOND_BARRIER 6
#define ADDREM_STAGPHASES_FIRST_BARRIER 7
#define ADDREM_STAGPHASES_SECOND_BARRIER 8
#define DOUBLE_REDUCE_FIRST_BARRIER 9
#define DOUBLE_REDUCE_SECOND_BARRIER 10
#define FLOAT_128_REDUCE_FIRST_BARRIER 11
#define FLOAT_128_REDUCE_SECOND_BARRIER 12
#define INTERNAL_NISSA_MALLOC_FIRST_BARRIER 13
#define INTERNAL_NISSA_MALLOC_SECOND_BARRIER 14
#define INTERNAL_NISSA_FREE_FIRST_BARRIER 15
#define INTERNAL_NISSA_FREE_SECOND_BARRIER 16
#define BUFFERED_COMM_WAIT_BARRIER 17
#define BUFFERED_COMM_LX_SENDING_BUF_FILL_BARRIER 18
#define BUFFERED_COMM_EV_OR_OD_SENDING_BUF_FILL_BARRIER 19
#define BUFFERED_COMM_EV_AND_OD_SENDING_BUF_FILL_BARRIER 20
#define INTERNAL_VECTOR_RESET_FIRST_BARRIER 21
#define INTERNAL_VECTOR_RESET_SECOND_BARRIER 22
#define HMC_SCALE_BARRIER 23
#define WILSON_STAPLE_BARRIER 24
#define CONTRACT_BARRIER 25
#define VECTOR_COPY_FIRST_BARRIER 26
#define VECTOR_COPY_SECOND_BARRIER 27
#define HOPPING_MATRIX_APPLICATION_BARRIER 28
#define REMAP_BARRIER 29
#define ATOMIC_FIRST_BARRIER 100
#define ATOMIC_SECOND_BARRIER 101

//////////////////////////////////////////////////////////////////////////////////////

#define GET_THREAD_ID() int thread_id=omp_get_thread_num()

#define THREAD_BARRIER_FORCE(a) thread_barrier(a,true)

#define IS_MASTER_THREAD (!thread_id)

#define NISSA_PARALLEL_LOOP(INDEX,EXT_START,EXT_END)			\
  for(int WORKLOAD=EXT_END-EXT_START,					\
	NCHUNKS=(thread_pool_locked)?1:nthreads,			\
	THREAD_LOAD=(WORKLOAD+NCHUNKS-1)/NCHUNKS,			\
	START=EXT_START+thread_id*THREAD_LOAD,				\
	END=START+THREAD_LOAD< EXT_END ? START+THREAD_LOAD : EXT_END,	\
	INDEX=START;INDEX<END;INDEX++)

//////////////////////////////////////////////////////////////////////////////////////

#define THREAD_ATOMIC_EXEC(inst)					\
  {thread_barrier(ATOMIC_FIRST_BARRIER);inst;thread_barrier(ATOMIC_SECOND_BARRIER);}

//////////////////////////////////////////////////////////////////////////////////////

//external argument to exchange info between function and worker
#define EXTERNAL_ARG(FUNC_NAME,LINE,ARG_TYPE,ARG) ARG_TYPE NAME4(FUNC_NAME,LINE,EXT_ARG,ARG);
#define EXPORT_ARG(FUNC_NAME,LINE,ARG) NAME4(FUNC_NAME,LINE,EXT_ARG,ARG)=ARG;
#define IMPORT_ARG(FUNC_NAME,LINE,ARG_TYPE,ARG) ARG_TYPE ARG=NAME4(FUNC_NAME,LINE,EXT_ARG,ARG);

//////////////////////////////////////////////////////////////////////////////////////

//headers: external parameters
#define THREADABLE_FUNCTION_0ARG_EXTERNAL_ARGS(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_1ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1)	\
  EXTERNAL_ARG(FUNC_NAME,LINE,AT1,A1)					\
  THREADABLE_FUNCTION_0ARG_EXTERNAL_ARGS(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_2ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2) \
  EXTERNAL_ARG(FUNC_NAME,LINE,AT2,A2)					\
  THREADABLE_FUNCTION_1ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1)
#define THREADABLE_FUNCTION_3ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3) \
  EXTERNAL_ARG(FUNC_NAME,LINE,AT3,A3)					\
  THREADABLE_FUNCTION_2ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2)
#define THREADABLE_FUNCTION_4ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4) \
  EXTERNAL_ARG(FUNC_NAME,LINE,AT4,A4)					\
  THREADABLE_FUNCTION_3ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3)
#define THREADABLE_FUNCTION_5ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) \
  EXTERNAL_ARG(FUNC_NAME,LINE,AT5,A5)					\
  THREADABLE_FUNCTION_4ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4)
#define THREADABLE_FUNCTION_6ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) \
  EXTERNAL_ARG(FUNC_NAME,LINE,AT6,A6)					\
  THREADABLE_FUNCTION_5ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5)
#define THREADABLE_FUNCTION_7ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) \
  EXTERNAL_ARG(FUNC_NAME,LINE,AT7,A7)					\
  THREADABLE_FUNCTION_6ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6)
#define THREADABLE_FUNCTION_8ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) \
  EXTERNAL_ARG(FUNC_NAME,LINE,AT8,A8)					\
  THREADABLE_FUNCTION_7ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7)
#define THREADABLE_FUNCTION_9ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) \
  EXTERNAL_ARG(FUNC_NAME,LINE,AT9,A9)					\
  THREADABLE_FUNCTION_8ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8)
#define THREADABLE_FUNCTION_10ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) \
  EXTERNAL_ARG(FUNC_NAME,LINE,AT10,A10)					\
  THREADABLE_FUNCTION_9ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9)
#define THREADABLE_FUNCTION_11ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11) \
  EXTERNAL_ARG(FUNC_NAME,LINE,AT11,A11)					\
  THREADABLE_FUNCTION_10ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10)

///////////////////////////////////////////////////////////////////////////////////////

//external function: exportation (last line is most external)
#define THREADABLE_FUNCTION_0ARG_EXPORT(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_1ARG_EXPORT(FUNC_NAME,LINE,AT1,A1)		\
  THREADABLE_FUNCTION_0ARG_EXPORT(FUNC_NAME,LINE)				\
  EXPORT_ARG(FUNC_NAME,LINE,A1)
#define THREADABLE_FUNCTION_2ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2)	\
  THREADABLE_FUNCTION_1ARG_EXPORT(FUNC_NAME,LINE,AT1,A1)			\
  EXPORT_ARG(FUNC_NAME,LINE,A2)
#define THREADABLE_FUNCTION_3ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3)	\
  THREADABLE_FUNCTION_2ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2)		\
  EXPORT_ARG(FUNC_NAME,LINE,A3)
#define THREADABLE_FUNCTION_4ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4) \
  THREADABLE_FUNCTION_3ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3)	\
  EXPORT_ARG(FUNC_NAME,LINE,A4)
#define THREADABLE_FUNCTION_5ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) \
  THREADABLE_FUNCTION_4ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4) \
  EXPORT_ARG(FUNC_NAME,LINE,A5)
#define THREADABLE_FUNCTION_6ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) \
  THREADABLE_FUNCTION_5ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) \
  EXPORT_ARG(FUNC_NAME,LINE,A6)
#define THREADABLE_FUNCTION_7ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) \
  THREADABLE_FUNCTION_6ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) \
  EXPORT_ARG(FUNC_NAME,LINE,A7)
#define THREADABLE_FUNCTION_8ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) \
  THREADABLE_FUNCTION_7ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) \
  EXPORT_ARG(FUNC_NAME,LINE,A8)
#define THREADABLE_FUNCTION_9ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) \
  THREADABLE_FUNCTION_8ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) \
  EXPORT_ARG(FUNC_NAME,LINE,A9)
#define THREADABLE_FUNCTION_10ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) \
  THREADABLE_FUNCTION_9ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) \
  EXPORT_ARG(FUNC_NAME,LINE,A10)
#define THREADABLE_FUNCTION_11ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11) \
  THREADABLE_FUNCTION_10ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) \
  EXPORT_ARG(FUNC_NAME,LINE,A11)

/////////////////////////////////////////////////////////////////////////////////////////

//body: issue worker and reimport (last line is again the most external, so this is why we have to split)
#define THREADABLE_FUNCTION_TRAMPOLINE_PT1()				\
  if(nthreads>1 && thread_pool_locked)					\
    {
#define THREADABLE_FUNCTION_TRAMPOLINE_PT2(FUNC_NAME,LINE)		\
      start_threaded_function(NAME3(FUNC_NAME,LINE,SUMMONER),#FUNC_NAME);	\
    }									\
  else

#define THREADABLE_FUNCTION_0ARG_TRAMPOLINE(FUNC_NAME,LINE)	\
  THREADABLE_FUNCTION_TRAMPOLINE_PT1()				\
  THREADABLE_FUNCTION_0ARG_EXPORT(FUNC_NAME,LINE)			\
  THREADABLE_FUNCTION_TRAMPOLINE_PT2(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_1ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1)	\
  THREADABLE_FUNCTION_TRAMPOLINE_PT1()					\
  THREADABLE_FUNCTION_1ARG_EXPORT(FUNC_NAME,LINE,AT1,A1)			\
  THREADABLE_FUNCTION_TRAMPOLINE_PT2(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_2ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT1()					\
  THREADABLE_FUNCTION_2ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2)		\
  THREADABLE_FUNCTION_TRAMPOLINE_PT2(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_3ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT1()					\
  THREADABLE_FUNCTION_3ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3)	\
  THREADABLE_FUNCTION_TRAMPOLINE_PT2(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_4ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT1()					\
  THREADABLE_FUNCTION_4ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT2(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_5ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT1()					\
  THREADABLE_FUNCTION_5ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT2(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_6ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT1()					\
  THREADABLE_FUNCTION_6ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT2(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_7ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT1()					\
  THREADABLE_FUNCTION_7ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT2(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_8ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT1()					\
  THREADABLE_FUNCTION_8ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT2(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_9ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT1()					\
  THREADABLE_FUNCTION_9ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT2(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_10ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT1()					\
  THREADABLE_FUNCTION_10ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT2(FUNC_NAME,LINE)
#define THREADABLE_FUNCTION_11ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT1()					\
  THREADABLE_FUNCTION_11ARG_EXPORT(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11) \
  THREADABLE_FUNCTION_TRAMPOLINE_PT2(FUNC_NAME,LINE)

/////////////////////////////////////////// summoner ///////////////////////////////////////////////////

#define THREADABLE_FUNCTION_0ARG_SUMMONER(FUNC_NAME,LINE)	\
  void NAME3(FUNC_NAME,LINE,SUMMONER)()				\
  {								\
    FUNC_NAME();						\
  }
#define THREADABLE_FUNCTION_1ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1)	\
  void NAME3(FUNC_NAME,LINE,SUMMONER)()					\
  {									\
    IMPORT_ARG(FUNC_NAME,LINE,AT1,A1);					\
    FUNC_NAME(A1);							\
  }
#define THREADABLE_FUNCTION_2ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2)	\
  void NAME3(FUNC_NAME,LINE,SUMMONER)()					\
  {									\
    IMPORT_ARG(FUNC_NAME,LINE,AT1,A1);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT2,A2);					\
    FUNC_NAME(A1,A2);							\
  }
#define THREADABLE_FUNCTION_3ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3) \
  void NAME3(FUNC_NAME,LINE,SUMMONER)()					\
  {									\
    IMPORT_ARG(FUNC_NAME,LINE,AT1,A1);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT2,A2);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT3,A3);					\
    FUNC_NAME(A1,A2,A3);							\
  }
#define THREADABLE_FUNCTION_4ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4) \
  void NAME3(FUNC_NAME,LINE,SUMMONER)()					\
  {									\
    IMPORT_ARG(FUNC_NAME,LINE,AT1,A1);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT2,A2);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT3,A3);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT4,A4);					\
    FUNC_NAME(A1,A2,A3,A4);						\
  }
#define THREADABLE_FUNCTION_5ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) \
  void NAME3(FUNC_NAME,LINE,SUMMONER)()					\
  {									\
    IMPORT_ARG(FUNC_NAME,LINE,AT1,A1);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT2,A2);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT3,A3);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT4,A4);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT5,A5);					\
    FUNC_NAME(A1,A2,A3,A4,A5);						\
  }
#define THREADABLE_FUNCTION_6ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) \
  void NAME3(FUNC_NAME,LINE,SUMMONER)()					\
  {									\
    IMPORT_ARG(FUNC_NAME,LINE,AT1,A1);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT2,A2);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT3,A3);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT4,A4);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT5,A5);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT6,A6);					\
    FUNC_NAME(A1,A2,A3,A4,A5,A6);					\
  }
#define THREADABLE_FUNCTION_7ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) \
  void NAME3(FUNC_NAME,LINE,SUMMONER)()					\
  {									\
    IMPORT_ARG(FUNC_NAME,LINE,AT1,A1);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT2,A2);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT3,A3);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT4,A4);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT5,A5);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT6,A6);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT7,A7);					\
    FUNC_NAME(A1,A2,A3,A4,A5,A6,A7);					\
  }
#define THREADABLE_FUNCTION_8ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) \
  void NAME3(FUNC_NAME,LINE,SUMMONER)()					\
  {									\
    IMPORT_ARG(FUNC_NAME,LINE,AT1,A1);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT2,A2);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT3,A3);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT4,A4);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT5,A5);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT6,A6);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT7,A7);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT8,A8);					\
    FUNC_NAME(A1,A2,A3,A4,A5,A6,A7,A8);					\
  }
#define THREADABLE_FUNCTION_9ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) \
  void NAME3(FUNC_NAME,LINE,SUMMONER)()					\
  {									\
    IMPORT_ARG(FUNC_NAME,LINE,AT1,A1);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT2,A2);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT3,A3);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT4,A4);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT5,A5);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT6,A6);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT7,A7);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT8,A8);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT9,A9);					\
    FUNC_NAME(A1,A2,A3,A4,A5,A6,A7,A8,A9);				\
  }
#define THREADABLE_FUNCTION_10ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) \
  void NAME3(FUNC_NAME,LINE,SUMMONER)()					\
  {									\
    IMPORT_ARG(FUNC_NAME,LINE,AT1,A1);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT2,A2);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT3,A3);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT4,A4);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT5,A5);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT6,A6);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT7,A7);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT8,A8);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT9,A9);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT10,A10);				\
    FUNC_NAME(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10);				\
  }
#define THREADABLE_FUNCTION_11ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11) \
  void NAME3(FUNC_NAME,LINE,SUMMONER)()					\
  {									\
    IMPORT_ARG(FUNC_NAME,LINE,AT1,A1);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT2,A2);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT3,A3);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT4,A4);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT5,A5);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT6,A6);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT7,A7);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT8,A8);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT9,A9);					\
    IMPORT_ARG(FUNC_NAME,LINE,AT10,A10);				\
    IMPORT_ARG(FUNC_NAME,LINE,AT11,A11);				\
    FUNC_NAME(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11);			\
  }

//////////////////////////////////////////////////////////////////////////////////////////

//threadable function with 0 arguments
#define THREADABLE_FUNCTION_0ARG_INSIDE(FUNC_NAME,LINE)			\
  THREADABLE_FUNCTION_0ARG_EXTERNAL_ARGS(FUNC_NAME,LINE)		\
  void FUNC_NAME();							\
  THREADABLE_FUNCTION_0ARG_SUMMONER(FUNC_NAME,LINE)			\
  void FUNC_NAME(){							\
    THREADABLE_FUNCTION_0ARG_TRAMPOLINE(FUNC_NAME,LINE)

//threadable function with 1 arguments
#define THREADABLE_FUNCTION_1ARG_INSIDE(FUNC_NAME,LINE,AT1,A1)		\
  THREADABLE_FUNCTION_1ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1)		\
  void FUNC_NAME(AT1 A1);						\
  THREADABLE_FUNCTION_1ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1)		\
  void FUNC_NAME(AT1 A1){						\
    THREADABLE_FUNCTION_1ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1)

//threadable function with 2 arguments
#define THREADABLE_FUNCTION_2ARG_INSIDE(FUNC_NAME,LINE,AT1,A1,AT2,A2)	\
  THREADABLE_FUNCTION_2ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2)	\
  void FUNC_NAME(AT1 A1,AT2 A2);					\
  THREADABLE_FUNCTION_2ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2)	\
  void FUNC_NAME(AT1 A1,AT2 A2){					\
    THREADABLE_FUNCTION_2ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2)

//threadable function with 3 arguments
#define THREADABLE_FUNCTION_3ARG_INSIDE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3) \
  THREADABLE_FUNCTION_3ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3);					\
  THREADABLE_FUNCTION_3ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3){					\
    THREADABLE_FUNCTION_3ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3)

//threadable function with 4 arguments
#define THREADABLE_FUNCTION_4ARG_INSIDE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4) \
  THREADABLE_FUNCTION_4ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4);				\
  THREADABLE_FUNCTION_4ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4){				\
    THREADABLE_FUNCTION_4ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4)

//threadable function with 5 arguments
#define THREADABLE_FUNCTION_5ARG_INSIDE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) \
  THREADABLE_FUNCTION_5ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5);			\
  THREADABLE_FUNCTION_5ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5){			\
    THREADABLE_FUNCTION_5ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5)

//threadable function with 6 arguments
#define THREADABLE_FUNCTION_6ARG_INSIDE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) \
  THREADABLE_FUNCTION_6ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6);		\
  THREADABLE_FUNCTION_6ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6){		\
    THREADABLE_FUNCTION_6ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6)

//threadable function with 7 arguments
#define THREADABLE_FUNCTION_7ARG_INSIDE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) \
  THREADABLE_FUNCTION_7ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7);	\
  THREADABLE_FUNCTION_7ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7){	\
    THREADABLE_FUNCTION_7ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7)
  
  //threadable function with 8 arguments
#define THREADABLE_FUNCTION_8ARG_INSIDE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) \
  THREADABLE_FUNCTION_8ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8); \
  THREADABLE_FUNCTION_8ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8){ \
    THREADABLE_FUNCTION_8ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8)

//threadable function with 9 arguments
#define THREADABLE_FUNCTION_9ARG_INSIDE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) \
  THREADABLE_FUNCTION_9ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8,AT9 A9); \
  THREADABLE_FUNCTION_9ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8,AT9 A9){ \
    THREADABLE_FUNCTION_9ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9)

//threadable function with 10 arguments
#define THREADABLE_FUNCTION_10ARG_INSIDE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) \
  THREADABLE_FUNCTION_10ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8,AT9 A9,AT10 A10); \
  THREADABLE_FUNCTION_10ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8,AT9 A9,AT10 A10){ \
    THREADABLE_FUNCTION_10ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10)

//threadable function with 11 arguments
#define THREADABLE_FUNCTION_11ARG_INSIDE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11) \
  THREADABLE_FUNCTION_11ARG_EXTERNAL_ARGS(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8,AT9 A9,AT10 A10,AT11 A11); \
  THREADABLE_FUNCTION_11ARG_SUMMONER(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11) \
  void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8,AT9 A9,AT10 A10,AT11 A11){ \
    THREADABLE_FUNCTION_11ARG_TRAMPOLINE(FUNC_NAME,LINE,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11)

///////////////////////////////////////////////////////////////////////////////////////////////////

#define THREADABLE_FUNCTION_0ARG(FUNC_NAME)		\
  THREADABLE_FUNCTION_0ARG_INSIDE(FUNC_NAME,__LINE__)
#define THREADABLE_FUNCTION_1ARG(FUNC_NAME,AT1,A1)		\
  THREADABLE_FUNCTION_1ARG_INSIDE(FUNC_NAME,__LINE__,AT1,A1)
#define THREADABLE_FUNCTION_2ARG(FUNC_NAME,AT1,A1,AT2,A2)		\
  THREADABLE_FUNCTION_2ARG_INSIDE(FUNC_NAME,__LINE__,AT1,A1,AT2,A2)
#define THREADABLE_FUNCTION_3ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3)	\
  THREADABLE_FUNCTION_3ARG_INSIDE(FUNC_NAME,__LINE__,AT1,A1,AT2,A2,AT3,A3)
#define THREADABLE_FUNCTION_4ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4)	\
  THREADABLE_FUNCTION_4ARG_INSIDE(FUNC_NAME,__LINE__,AT1,A1,AT2,A2,AT3,A3,AT4,A4)
#define THREADABLE_FUNCTION_5ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) \
  THREADABLE_FUNCTION_5ARG_INSIDE(FUNC_NAME,__LINE__,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5)
#define THREADABLE_FUNCTION_6ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) \
  THREADABLE_FUNCTION_6ARG_INSIDE(FUNC_NAME,__LINE__,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6)
#define THREADABLE_FUNCTION_7ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) \
  THREADABLE_FUNCTION_7ARG_INSIDE(FUNC_NAME,__LINE__,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7)
#define THREADABLE_FUNCTION_8ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) \
  THREADABLE_FUNCTION_8ARG_INSIDE(FUNC_NAME,__LINE__,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8)
#define THREADABLE_FUNCTION_9ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) \
  THREADABLE_FUNCTION_9ARG_INSIDE(FUNC_NAME,__LINE__,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9)
#define THREADABLE_FUNCTION_10ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) \
  THREADABLE_FUNCTION_10ARG_INSIDE(FUNC_NAME,__LINE__,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10)
#define THREADABLE_FUNCTION_11ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11) \
  THREADABLE_FUNCTION_11ARG_INSIDE(FUNC_NAME,__LINE__,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11)

#endif
