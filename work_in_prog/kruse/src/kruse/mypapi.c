#include "mypapi.h"
#include "config.h"


#if PAPI

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
//#include <common/alignment.h>
//#include <papi.h> //Provides definitions for base PAPI function
//#include <spi/bgp_SPI.h> //The Blue Gene/P substrate for PAPI is based on the UPC interfaces documented in this header file.
//#include <papiStdEventDefs.h> //Provides a listing of all standard PAPI presets. Please note that not all of these are available on the Blue Gene/P platform
//#include <linux-bgp-native-events.h> //Provides a listing of all available counters native to Blue Gene.
#include <omp.h>

//#include <linux-bgq.h>
//#include <fpapi.h>
#include <upci/events.h>
#include <mpi.h>
#include <assert.h>
#include <kernel/location.h>
#include <bgpm/include/bgpm.h>

//void FPUArith(void); //This method does various calculations which should saturate many of the counters
//void List_PAPI_Events(const int pEventSet, int* pEvents, int* xNumEvents);
//void Print_Native_Counters();
//void Print_Native_Counters_via_Buffer(const BGP_UPC_Read_Counters_Struct_t* pBuffer);
//void Print_Native_Counters_for_PAPI_Counters(const int pEventSet);
//void Print_Native_Counters_for_PAPI_Counters_From_List(const int* pEvents, const int pNumEvents);
//void Print_PAPI_Counters(const int pEventSet, const long long* pCounters);
//void Print_PAPI_Counters_From_List(const int* pEventList, const int pNumEvents, const long long* pCounters);
void Print_Counters(const int pEventSet);
//void Print_PAPI_Events(const int pEventSet);
long long getMyPapiValue(const int eventNum);

#define lengthof(X) (sizeof(X)/sizeof((X)[0]))

int PAPI_Events[256];
long long PAPI_Counters[256];
//int xEventSet=PAPI_NULL;

long long xCyc;
long long xNsec;
double xNow;
static double xWtime;
static double xOmpTime;

extern int g_proc_id;

static double now2(){
   struct timeval t; double f_t;
   gettimeofday(&t, NULL);
   f_t = t.tv_usec; f_t = f_t/1000000.0; f_t +=t.tv_sec;
   return f_t;
}

#define PAPI_ERROR(cmd)                                                                                      \
	do {                                                                                                       \
		int RC = (cmd);                                                                                       \
		if (RC != PAPI_OK) {                                                                                   \
			 fprintf(stderr, "MK_PAPI call failed with code %d at line %d on MPI rank %d thread %d: %s\n", RC, __LINE__,  g_proc_id, Kernel_ProcessorID(), TOSTRING(cmd)); \
		}                                                                                                        \
	} while (0)


#define BGPM_ERROR(cmd)                                                                                      \
	do {                                                                                                       \
		int RC = (cmd);                                                                                       \
		if (RC) {                                                                                   \
			 fprintf(stderr, "MK_BGPM call failed with code %d at line %d on MPI rank %d thread %d: %s\n", RC, __LINE__,  g_proc_id, Kernel_ProcessorID(), TOSTRING(cmd)); \
		}                                                                                                        \
	} while (0)


//static int PAPI_add_native_event(int EventSet, int EventCode) {
	// For some strange reason (i.e. unknown to me), the numbers from events.h are off by one
//	return PAPI_add_event(EventSet, PAPI_NATIVE_MASK | (EventCode-1));
//}


mypapi_counters mypapi_merge_counters(mypapi_counters *counters1, mypapi_counters *counters2) {
	if (!counters1->init && !counters2->init) {
		mypapi_counters emptyresult = { 0 };
		return emptyresult;
	}

	mypapi_counters result;
	if (!counters2->init) {
		//if (g_proc_id == 0 && omp_get_thread_num() == 0)
		//		fprintf(stderr, "MK Take counter1: %llu\n", counters1->native[PEVT_CYCLES]);
		return *counters1;
	} else if (!counters1->init) {
		//if (g_proc_id == 0 && omp_get_thread_num() == 0)
		//		fprintf(stderr, "MK Take counter2: %llu\n", counters2->native[PEVT_CYCLES]);
		return *counters2;
	} else {
		if (counters1->set == counters2->set)
			result.set = counters1->set;
		else
			result.set = -1;

		if (counters1->eventset == counters2->eventset)
			result.eventset = counters1->eventset;
		else
			result.eventset = -1;

		if (counters1->threadid == counters2->threadid)
			result.threadid = counters1->threadid;
		else
			result.threadid = -1;

		if (counters1->coreid == counters2->coreid)
			result.coreid = counters1->coreid;
		else
			result.coreid = -1;

		if (counters1->smtid == counters2->smtid)
			result.smtid = counters1->smtid;
		else
			result.smtid = -1;

		if (counters1->ompid == counters2->ompid)
			result.ompid = counters1->ompid;
		else
			result.ompid = -1;

		//for (int i = 0; i < lengthof(counters1->preset); i += 1) {
		//	result.preset[i] = counters1->preset[i] + counters2->preset[i];
		//}
		for (int i = 0; i < lengthof(counters1->native); i += 1) {
			uint64_t c1 = counters1->native[i];
			uint64_t c2 = counters2->native[i];
			uint64_t merged;
			switch (i) {
			case PEVT_CYCLES:
				// Count PEVT_CYCLES just once
				if (counters1->set > 0)
					c1 = 0;
				if (counters2->set > 0)
					c2 = 0;
				// Fallthrough
			default:
				merged = c1 + c2;
				break;
			}
			result.native[i] = merged;
			result.active[i] = counters1->active[i] || counters2->active[i];
		}

		result.corecycles = ((counters1->set > 0) ?  0 : counters1->corecycles) + ((counters2->set > 0) ? 0 : counters2->corecycles);
		result.nodecycles = ((counters1->set > 0) ?  0 : counters1->nodecycles) + ((counters2->set > 0) ? 0 : counters2->nodecycles);

		if (g_proc_id == 0 && omp_get_thread_num() == 0) {
			//fprintf(stderr, "MK Merge result: %llu\n", result.native[PEVT_CYCLES]);
		}
		result.init = true;
	}

	assert(result.init);
	return result;
}


static int mypapi_eventsets[64][MYPAPI_SETS];

static unsigned long int mypapi_getthreadid() {
	return Kernel_ProcessorID();
}


static int PuEventSets[MYPAPI_SETS][64] = {0};
static int L2EventSet[MYPAPI_SETS] = {0};


static mypapi_counters mypapi_bgpm_read(int eventset, int set) {
	mypapi_counters result = { 0 };
	result.set = set;
	result.eventset = eventset;
	result.threadid = Kernel_ProcessorID(); /* 0..63, Kernel_ProcessorCoreID() << 2 + Kernel_ProcessorThreadID() */
	result.coreid = Kernel_ProcessorCoreID(); /* 0..15 */
	result.smtid = Kernel_ProcessorThreadID(); /* 0..3 */
	result.ompid = omp_get_thread_num();

	int numEvts = Bgpm_NumEvents(eventset);
	assert(numEvts > 0);

	uint64_t cnt;
	for (int i = 0; i < numEvts; i += 1) {
		int eventid = Bgpm_GetEventId(eventset, i);

		switch (eventid) {
		case PEVT_L1P_BAS_LU_STALL_SRT:
		case PEVT_L1P_BAS_LU_STALL_SRT_CYC:
		case PEVT_L1P_BAS_LU_STALL_MMIO_DCR:
		case PEVT_L1P_BAS_LU_STALL_MMIO_DCR_CYC:
		case PEVT_L1P_BAS_LU_STALL_STRM_DET:
		case PEVT_L1P_BAS_LU_STALL_STRM_DET_CYC:
		case PEVT_L1P_BAS_LU_STALL_LIST_RD:
		case PEVT_L1P_BAS_LU_STALL_LIST_RD_CYC:
		case PEVT_L1P_BAS_ST:
		case PEVT_L1P_BAS_LU_STALL_LIST_WRT:
		case PEVT_L1P_BAS_LU_STALL_LIST_WRT_CYC:

		case PEVT_L1P_STRM_LINE_ESTB:
		case PEVT_L1P_STRM_HIT_FWD:
		case PEVT_L1P_STRM_L1_HIT_FWD:
		case PEVT_L1P_STRM_EVICT_UNUSED:
		case PEVT_L1P_STRM_EVICT_PART_USED:
		case PEVT_L1P_STRM_REMOTE_INVAL_MATCH:
		case PEVT_L1P_STRM_DONT_CACHE:
		case PEVT_L1P_STRM_LINE_ESTB_ALL_LIST:
			// Per core counters
			if (result.smtid != 0)
				continue;
			break;
		case PEVT_L2_HITS:
		case PEVT_L2_MISSES:
		case PEVT_L2_FETCH_LINE:
		case PEVT_L2_STORE_LINE:
		case PEVT_L2_PREFETCH:
		case PEVT_L2_STORE_PARTIAL_LINE:
			// Per node counters, store just one event
			if (result.threadid != 0)
				continue;
			break;
		}

		uint64_t cnt;
		BGPM_ERROR(Bgpm_ReadEvent(eventset, i, &cnt));
		result.native[eventid] = cnt;
		result.active[eventid] = true;
	}

	// Cycles are the some on the whole node, but shared between smt threads
	// Just one of 2 threads on the same core is issued, so we count cycles one per core
	if (result.smtid == 0)
		result.corecycles = result.native[PEVT_CYCLES];
	else
		result.corecycles = 0;

	if (result.threadid == 0)
		result.nodecycles = result.native[PEVT_CYCLES];
	else
		result.nodecycles = 0;

	result.init = true;
	if (g_proc_id == 0 && omp_get_thread_num() == 0) {
		//mypapi_print_counters(&result);
	}
	return result;
}


void mypapi_init() {
	//assert(omp_get_thread_num() == 0);
	//if ((g_proc_id == 0) && (omp_get_thread_num() == 0))
	//	fprintf(stderr, "MK_Init mypapi\n");

//#pragma omp parallel
	{
		BGPM_ERROR(Bgpm_Init(BGPM_MODE_SWDISTRIB));

		int tid = Kernel_ProcessorID();
		int cid = Kernel_ProcessorCoreID();
		int sid = Kernel_ProcessorThreadID();
		for (int i = 0; i < MYPAPI_SETS; i += 1) {
			PuEventSets[i][tid] = Bgpm_CreateEventSet();
			assert(PuEventSets[i][tid] >= 0);
		}

		if (tid == 0) {
			for (int i = 0; i < MYPAPI_SETS; i += 1) {
				L2EventSet[i] = Bgpm_CreateEventSet();
				assert(L2EventSet[i] >= 0);
			}
		}

		int j = 0;
		{
			int pues = PuEventSets[j][tid];
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_CYCLES));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_INST_ALL));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_IU_IL1_MISS_CYC));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_IU_IBUFF_EMPTY_CYC));
			//BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_LIST_MISMATCH)); // Conflict set #3 // core address does not match a list address

			if (tid == 0) {
				int l2es = L2EventSet[j];
				BGPM_ERROR(Bgpm_AddEvent(l2es, PEVT_L2_HITS));
				BGPM_ERROR(Bgpm_AddEvent(l2es, PEVT_L2_MISSES));
				BGPM_ERROR(Bgpm_AddEvent(l2es, PEVT_L2_FETCH_LINE)); // L2 lines loaded from main memory
				BGPM_ERROR(Bgpm_AddEvent(l2es, PEVT_L2_STORE_LINE)); // L2 lines stored to main memory
				BGPM_ERROR(Bgpm_AddEvent(l2es, PEVT_L2_PREFETCH));
				BGPM_ERROR(Bgpm_AddEvent(l2es, PEVT_L2_STORE_PARTIAL_LINE));
			}
		}

		j += 1;
		{
			int pues = PuEventSets[j][tid];
			//BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_CYCLES));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_LSU_COMMIT_STS)); // Number of completed store commands.
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_LSU_COMMIT_LD_MISSES));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_BAS_STRM_LINE_ESTB)); // Conflict set #1 // Lines established for stream prefetch
			//BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_LIST_SKIP)); // Conflict set #2 // core address matched a non head of queue list address (per thread)
			//BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_LIST_STARTED)); // Conflict set #4 // List prefetch process was started
			//BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_LIST_OVF_MEM)); // Conflict set #5 // Written pattern exceeded allocated buffer

			if (tid == 0) {
				int l2es = L2EventSet[j];
				BGPM_ERROR(Bgpm_DeleteEventSet(l2es));
				L2EventSet[j] = -1;
			}
		}

		j += 1;
		{
			int pues = PuEventSets[j][tid];
			//BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_CYCLES));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_LSU_COMMIT_DCBT_MISSES)); // Number of completed dcbt[st][ls][ep] commands that missed the L1 Data Cache.
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_LSU_COMMIT_DCBT_HITS)); // Number of completed dcbt[st][ls][ep] commands that missed the L1 Data Cache.
			//BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_LIST_CMP)); // Conflict set #1 // core address was compared against list
			//BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_LIST_ABANDON)); // Conflict set #5 // A2 loads mismatching pattern resulted in abandoned list prefetch
			//BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_LIST_CMP_OVRUN_PREFCH)); // Conflict set #6 // core address advances faster than prefetch lines can be established dropping prefetches

			if (tid == 0) {
				int l2es = L2EventSet[j];
				BGPM_ERROR(Bgpm_DeleteEventSet(l2es));
				L2EventSet[j] = -1;
			}
		}

		j += 1;
		{
			int pues = PuEventSets[j][tid];
			//BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_CYCLES));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_LSU_COMMIT_CACHEABLE_LDS)); // Number of completed cache-able load commands. (without dcbt)
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_INST_XU_ALL)); // All XU instructions completed (instructions which use A2 FX unit - UPC_P_XU_OGRP_*).
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_INST_QFPU_ALL)); // Count all completed instructions which processed by the QFPU unit (UPC_P_AXU_OGRP_*)
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_BAS_MISS));// Conflict set #4
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_BAS_HIT)); // Conflict set #2 // Hits in prefetch directory
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_BAS_LU_STALL_LIST_RD_CYC)); // Conflict Set #1 // Cycles lookup was held while list fetched addresses

			if (tid == 0) {
				int l2es = L2EventSet[j];
				BGPM_ERROR(Bgpm_DeleteEventSet(l2es));
				L2EventSet[j] = -1;
			}
		}

		j += 1;
		{
			int pues = PuEventSets[j][tid];
			//BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_CYCLES));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_IU_IS1_STALL_CYC)); // Register Dependency Stall
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_IU_IS2_STALL_CYC)); // Instruction Issue Stall
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_STRM_EVICT_UNUSED)); // Conflict set #4 // per core // Lines fetched and never hit evicted
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_STRM_EVICT_PART_USED)); // Conflict set #5 // per core // Line fetched and only partially used evicted
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_STRM_LINE_ESTB)); // Conflict set #1 // per core // lines established for any reason and thread



			if (tid == 0) {
				int l2es = L2EventSet[j];
				BGPM_ERROR(Bgpm_DeleteEventSet(l2es));
				L2EventSet[j] = -1;
			}
		}

		j += 1;
		{
			int pues = PuEventSets[j][tid];
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_CYCLES));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_STRM_STRM_ESTB)); // Conflict set #1 // per thread // streams detected and established
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_STRM_HIT_LIST)); // Conflict set #4 // per thread // Hits for lines fetched by list engine

			if (tid == 0) {
				int l2es = L2EventSet[j];
				BGPM_ERROR(Bgpm_DeleteEventSet(l2es));
				L2EventSet[j] = -1;
			}
		}
	}
}


void mypapi_free() {
	BGPM_ERROR(Bgpm_Disable());
}


static double mypapi_wtime() {
	Personality_t personality;
	BGPM_ERROR(Kernel_GetPersonality(&personality, sizeof(Personality_t)));
	double freq = MEGA * personality.Kernel_Config.FreqMHz;
	long long cycles = GetTimeBase();
	return cycles / freq;
}

static int activeEventSet = -1;


static void mypapi_start_work(int i) {
	int tid = Kernel_ProcessorID();
	int cid = Kernel_ProcessorCoreID();
	int sid = Kernel_ProcessorThreadID();

	if (tid == 0) {
		xCyc = GetTimeBase();
		//xCyc = PAPI_get_real_cyc();
		xNsec = mypapi_wtime();
		//xNsec = PAPI_get_real_nsec();
		xNow = now2();
		xWtime = MPI_Wtime();
		xOmpTime = omp_get_wtime();

		activeEventSet = i;
	}

	int pues = PuEventSets[i][tid];
	BGPM_ERROR(Bgpm_Apply(pues));
	if (tid == 0) {
		int l2es = L2EventSet[i];
		if (l2es > 0) {
			BGPM_ERROR(Bgpm_Apply(l2es));
			BGPM_ERROR(Bgpm_Start(l2es));
		}
	}
	BGPM_ERROR(Bgpm_Start(pues));
}


void mypapi_start(int i) {
	if (omp_in_parallel()) {
		// Already in parallel, no need to start threads
		mypapi_start_work(i);
	} else {
		// Start counters in all threads
		#pragma omp parallel
		{
			mypapi_start_work(i);
		}
	}
}


void mypapi_print_counters(mypapi_counters *counters) {
	fprintf(stderr, "*******************************************************************************\n");
	fprintf(stderr, "Set=%d eventset=%d thread=%d core=%d smt=%d omp=%d\n", counters->set, counters->eventset, counters->threadid, counters->coreid, counters->smtid, counters->ompid);
	if (counters->corecycles) {
		fprintf(stderr, "%10llu = %-30s\n", counters->corecycles, "Core cycles");
	}
	for (int i = 0; i < lengthof(counters->native); i+=1) {
		if (counters->active[i]) {
			uint64_t val = counters->native[i];
			//const char *label = Bgpm_GetEventIdLabel(i);
			Bgpm_EventInfo_t info;
			BGPM_ERROR(Bgpm_GetEventIdInfo(i, &info));

			fprintf(stderr, "%10llu = %-30s (%s)\n", val, info.label, info.desc);
		}
	}
	fprintf(stderr, "*******************************************************************************\n");
}




void mypapi_stop_work(mypapi_counters *result) {
	int ompid = omp_get_thread_num();

	ompid = omp_get_thread_num();
	int tid = Kernel_ProcessorID();
	int cid = Kernel_ProcessorCoreID();
	int sid = Kernel_ProcessorThreadID();
	int i = activeEventSet;
	assert(i >= 0);
	assert(i < MYPAPI_SETS);

	int pues = PuEventSets[i][tid];
	BGPM_ERROR(Bgpm_Stop(pues));
	if (tid == 0) {
		int l2es = L2EventSet[i];
		if (l2es >= 0) {
			BGPM_ERROR(Bgpm_Stop(l2es));
		}
	}

	mypapi_counters local_result = mypapi_bgpm_read(pues, i);
	if (tid == 0) {
		int l2es = L2EventSet[i];
		if (l2es >= 0) {
			mypapi_counters local_result_l2 = mypapi_bgpm_read(l2es, i);
			local_result = mypapi_merge_counters(&local_result, &local_result_l2);
		}
	}

#pragma omp critical (mypapi)
	{
		*result = mypapi_merge_counters(result, &local_result);
	}
}


mypapi_counters mypapi_stop() {
	static mypapi_counters result; // static to force it shared between all threads, even if this func is called by all threads (i.e. in a #pragma omp parallel)
	if (omp_in_parallel()) {
		// Already in parallel, no need to start threads
		mypapi_stop_work(&result);
	} else {
		// Start counters in all threads
		#pragma omp parallel
		{
			mypapi_stop_work(&result);
		}
	}
	return result;
}

#if 0
/* List_PAPI_Events */
static void List_PAPI_Events(const int pEventSet, int* pEvents, int* pNumEvents) {
	int xRC = PAPI_list_events(pEventSet, pEvents, pNumEvents);
	if (xRC != PAPI_OK) {
		printf("FAILURE: PAPI_list_events failed, returned xRC=%d...\n", xRC);
		exit(1);
	}
	return;
}
#endif

#else

void mypapi_init(){}
void mypapi_start(int i){}
mypapi_counters mypapi_stop(){
	mypapi_counters dummy;
	dummy.init=false;
	return dummy;
}
void mypapi_free(){}

mypapi_counters mypapi_merge_counters(mypapi_counters *counters1, mypapi_counters *counters2) {
	mypapi_counters result = {0};
	return result;
}
void mypapi_print_counters(mypapi_counters *counters) {}



#endif

