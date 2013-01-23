/*
 * bgq_dispatch.c
 *
 *  Created on: Oct 15, 2012
 *      Author: meinersbur
 */

#define BGQ_DISPATCH_C_
#include "bgq_dispatch.h"

#include "bgq_qpx.h"

#if BGQ_QPX
#include <l2/barrier.h>
//#include <wu/wait.h>
#include <upci/upc_atomic.h>
//#include <hwi/include/bqc/A2_inlines.h>
//#include <time.h> // nanosleep() system call
#endif

#include "../global.h"

#include <omp.h>
#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

//static L2_Barrier_t barrier;

#ifndef NDEBUG
static char space[64*25+1]= {' '};
#endif

static volatile bgq_worker_func g_bgq_dispatch_func;
static void * volatile g_bgq_dispatch_arg;
static volatile bool g_bgq_dispatch_terminate;
static volatile bool g_bgq_dispatch_sync; // obsolete
//static volatile size_t g_bgq_dispatch_seq; // obsolete

static bool g_bgq_dispatch_pendingsync; // Accessed by master only


#if BGQ_QPX
static bool g_bgq_dispatch_barrier_initialized = false;
static L2_Barrier_t g_bgq_dispatch_barrier = L2_BARRIER_INITIALIZER;
#endif

static inline void bgq_thread_barrier() {
#if BGQ_QPX
	uint64_t savpri = Set_ThreadPriority_Low(); // Lower thread priority, so if busy waiting is used, do not impact other threads on core
	L2_Barrier(&g_bgq_dispatch_barrier, g_bgq_dispatch_threads);
	Restore_ThreadPriority(savpri);
#else
#pragma omp barrier
#endif
}


int bgq_parallel(bgq_master_func master_func, void *master_arg) {
#ifndef NDEBUG
	for (int i = 0; i < 64*25; i+=1)
		space[64*25] = ' ';
	space[64*25] = '\0';
#endif
	assert(!omp_in_parallel() && "This starts the parallel section, do not call it within one");
	g_bgq_dispatch_func = NULL;
	g_bgq_dispatch_arg = NULL;
	g_bgq_dispatch_terminate = false;
	g_bgq_dispatch_sync = false;
	//g_bgq_dispatch_seq = 0;
#if BGQ_QPX
	if (!g_bgq_dispatch_barrier_initialized) {
		Kernel_L2AtomicsAllocate(&g_bgq_dispatch_barrier, sizeof(g_bgq_dispatch_barrier));
		g_bgq_dispatch_barrier_initialized = true;
	}
#endif
	g_bgq_dispatch_threads = omp_get_max_threads();
#ifdef OMP
	omp_num_threads = 1/*omp_get_num_threads()*/; // For legacy linalg (it depends on whether nested parallelism is enabled)
#endif
	g_bgq_dispatch_inparallel = true;

	int master_result = 0;
	// We use OpenMP only to start the threads
	// Overhead of using OpenMP is too large
#pragma omp parallel
	{
		size_t tid = omp_get_thread_num(); // Or

		// Start workers
		if (tid != 0) {
			g_bgq_dispatch_pendingsync = true;
			bgq_worker();
		}

		// Start control program in master
		if (tid == 0) {
			master_result = master_func(master_arg);

			bgq_master_sync();
			// After program finishes, set flag so worker threads can terminate
			g_bgq_dispatch_func = NULL;
			g_bgq_dispatch_arg = NULL;
			g_bgq_dispatch_terminate = true;
			g_bgq_dispatch_sync = false;
#if BGQ_QPX
			mbar();
#else
#pragma omp flush
#endif

			// Wakeup workers to terminate
			bgq_worker();
		}

		//printf("%*sEXIT: tid=%u\n", (int)tid*25, "",(int)tid);
	}
	g_bgq_dispatch_inparallel = false;
	g_bgq_dispatch_threads = 0;
#ifdef OMP
	omp_num_threads = omp_get_max_threads();
#endif
	return master_result;
}


//static size_t count = 0;
//#pragma omp threadvar(count)

void bgq_worker() {
	//assert(omp_in_parallel() && "Call this inside #pragma omp parallel");
	assert(g_bgq_dispatch_threads == omp_get_num_threads());
	size_t threads = g_bgq_dispatch_threads;
	size_t tid = omp_get_thread_num(); // Or Kernel_ProcessorID()


	//assert((tid != 0) && "This function is for non-master threads only");
//size_t count = 0;
	while (true) {
		// Wait until every thread did its work
		// This doesn't need to be a barrier, waiting for submission of some work from the master is ok too
		// TODO: Hope OpenMP has a good implementation without busy waiting; if not, do some own work

		if (tid!=0) {
			// Guarantee that work is finished
			//TODO: can we implement this without this second barrier?
			// Required to ensure consistency of g_bgq_dispatch_sync, g_bgq_dispatch_terminate, g_bgq_dispatch_func
			bgq_thread_barrier(); // the sync barrier
		}


		// Master thread may write shared variables before this barrier
		bgq_thread_barrier();
#if BGQ_QPX
		mbar();
#else
#pragma omp flush
#endif
		// Worker threads read common variables after this barrier
		if (tid==0) {
			assert(!g_bgq_dispatch_pendingsync);
			g_bgq_dispatch_pendingsync = true;
		}

		//count += 1;
		// All threads should be here at the same time, including the master thread, which has issued some work, namely, calling a function

		if (g_bgq_dispatch_sync) {
			// This was meant just for synchronization between the threads, which already has been done
			//printf("%*sSYNC: tid=%u seq=%u\n", (int)tid*20, "", (int)tid, (int)g_bgq_dispatch_seq);
		} else if (g_bgq_dispatch_terminate) {
			// Exit program, or at least, the parallel section
			//printf("%*sTERM: tid=%u seq=%u\n", (int)tid*20, "",(int)tid, (int)g_bgq_dispatch_seq);
			return;
		} else {
			//printf("%*sCALL: tid=%u seq=%u\n", (int)tid*20, "",(int)tid, (int)g_bgq_dispatch_seq);
			assert(g_bgq_dispatch_func);
			void *arg = g_bgq_dispatch_arg;
			g_bgq_dispatch_func(arg, tid, threads); //TODO: Shuffle tid to load-balance work?
		}

		if (tid==0) {
			// Let master thread continue the program
			// Hint: master thread must call bgq_thread_barrier() sometime to release the workers from the following barrier
			return;
		}
		// All others, wait for the next command
	}
}


typedef struct {
	bgq_worker_func func;
	void *funcargs;
} bgq_master_adhoc_parms;

static int bgq_master_adhoc(void *arg_untyped) {
	bgq_master_adhoc_parms *args = arg_untyped;
	bgq_worker_func func = args->func;
	void *funcargs = args->funcargs;

	bgq_master_call(func, funcargs);
	return 0;
}

void bgq_master_call(bgq_worker_func func, void *arg) {
	assert(omp_get_thread_num()==0);
	assert(func);

	if (g_bgq_dispatch_inparallel) {
		bgq_master_sync();

		//g_bgq_dispatch_seq += 1;
		//printf("MASTER CALL seq=%d--------------------------------------------------------\n", (int)g_bgq_dispatch_seq);
		g_bgq_dispatch_func = func;
		g_bgq_dispatch_arg = arg;
		g_bgq_dispatch_terminate = false;
		g_bgq_dispatch_sync = false;
#if BGQ_QPX
	mbar();
#else
#pragma omp flush
#endif
		// Join the worker force
		bgq_worker();
	} else {
		// Not in a parallel section. Two possibilities:
		// 1. Execute sequentially
		// 2. Start and end a parallel section just for calling this func

		if (omp_in_parallel()) {
			// We are in another OpenMP construct, therefore call sequentially
			(*func)(arg, 0, 1);
		} else {
			bgq_master_adhoc_parms work = { func, arg };
			bgq_parallel(&bgq_master_adhoc, &work);
		}
	}
}


void bgq_master_sync() {
	assert(omp_get_thread_num()==0);

	if (g_bgq_dispatch_inparallel) {
		if (g_bgq_dispatch_pendingsync) {
			bgq_thread_barrier();
			//fflush(stdout);
			//printf("MASTER SYNC seq=%d--------------------------------------------------------\n", (int)g_bgq_dispatch_seq);
			//fflush(stdout);
			g_bgq_dispatch_pendingsync = false;
			return;
		} else {
			// Threads already at sync barrier
			return;
		}
	} else {
		// No sync necessary
	}
}


typedef struct {
	void *ptr;
	size_t size;
} bgq_memzero_work;

void bgq_memzero_worker(void *args_untyped, size_t tid, size_t threads) {
	bgq_memzero_work *args = args_untyped;
	char *ptr = args->ptr;
	size_t size = args->size;

	const size_t workload = size;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);

	char *beginLine = (ptr + begin);
	beginLine = (void*)((uintptr_t)beginLine & ~(BGQ_ALIGNMENT_L1-1));

	char *endLine  = (ptr + end);
	endLine = (void*)((uintptr_t)endLine & ~(BGQ_ALIGNMENT_L1-1));

	// Special cases
	//if (tid == 0)
	//	beginLine = ptr; /* rule redundant */
	if (tid == g_bgq_dispatch_threads-1)
		endLine = (ptr + size);

	if (beginLine < ptr)
		beginLine = ptr;
	if (beginLine >= endLine)
		return;

	assert(beginLine >= ptr);
	assert(endLine <= ptr + size);
	assert(endLine >= beginLine);

	size_t threadsize = (endLine-beginLine);
	memset(beginLine, 0x00, threadsize);
}

void bgq_master_memzero(void *ptr, size_t size) {
	if (size <= (1<<14)/*By empirical testing*/) {
		memset(ptr, 0x00, size);
	} else {
		bgq_master_sync();
		static bgq_memzero_work work;
		work.ptr = ptr;
		work.size = size;
		bgq_master_call(&bgq_memzero_worker, &work);
	}
}


typedef struct {
	void *ptrDst;
	void *ptrSrc;
	size_t size;
} bgq_memcopy_work;

void bgq_memcpy_worker(void *args_untyped, size_t tid, size_t threads) {
	bgq_memcopy_work *args = args_untyped;
	char *ptrDst = args->ptrDst;
	char *ptrSrc = args->ptrSrc;
	size_t size = args->size;

	const size_t workload = size;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);

	if (begin<end) {
		size_t count = end-begin;
		memcpy(ptrDst+begin, ptrSrc+begin, count);
	}
}

void bgq_master_memcpy(void *ptrDst, void *ptrSrc, size_t size) {
	if (size <= (1<<14)) {
		memcpy(ptrDst, ptrSrc, size);
	} else {
		bgq_master_sync();
		static bgq_memcopy_work work;
		work.ptrDst = ptrDst;
		work.ptrSrc = ptrSrc;
		work.size = size;
		bgq_master_call(&bgq_memcpy_worker, &work);
	}
}


typedef struct {
	bgq_mainlike_func func;
	int argc;
	char **argv;
} bgq_dispatch_mainlike_args;

int bgq_dispatch_callMainlike(void *args_untyped) {
	bgq_dispatch_mainlike_args *args = args_untyped;
	bgq_mainlike_func func = args->func;

	return (*func)(args->argc, args->argv);
}

int bgq_parallel_mainlike(bgq_mainlike_func func, int argc, char *argv[]) {
	bgq_dispatch_mainlike_args args = { func, argc, argv };
	return bgq_parallel(&bgq_dispatch_callMainlike, &args);
}

