#ifndef MYPAPI
#define MYPAPI


#include <stdbool.h>
#include <stdint.h>

#ifdef XLC
//#include <papi.h>
#include <upci/events.h>
#else
#define UPCI_NUM_EVENTS 5
#endif

#ifndef STRINGIFY
#define STRINGIFY(x) #x
#endif

#ifndef TOSTRING
#define TOSTRING(x) STRINGIFY(x)
#endif 

#define __Pragma(x) _Pragma(#x)

typedef struct {
	bool init;
	int set;
	int eventset;
	int threadid;
	int coreid;
	int smtid;
	int ompid;
	//long long preset[PAPI_END_idx];
	uint64_t native[UPCI_NUM_EVENTS];
	uint64_t corecycles;
	uint64_t nodecycles;
	bool active[UPCI_NUM_EVENTS];
	double secs;
} mypapi_counters;

typedef enum {
	pi_correct,
	pi_flopsref,
	pi_floppersite,
	pi_msecs,
	pi_cycpersite,

	pi_localrms,
	pi_globalrms,
	pi_avgovhtime,

	pi_cpi,
	pi_l1istalls,
	pi_axufraction,
	//pi_overhead,

	pi_is1stalls,
	pi_is2stalls,

	pi_hitinl1,
	pi_l1phitrate,
	pi_l2hitrate,
	pi_dcbthitrate,

	pi_detstreams,
	pi_l1pstreamunusedlines,

	pi_ramfetchrate,
	pi_ramstorerate,
	pi_ramstorepartial,
	pi_l2prefetch,

	__pi_COUNT,
	pi_corecpi,
	pi_instrpersite,
	pi_fxupersite,
	pi_flops,
	pi_l1pstreamhitinl1p,
	pi_hitinl1p,
	pi_l1pliststarted,
	pi_l1plistabandoned,
	pi_l1plistmismatch,
	pi_l1plistskips,
	pi_l1plistoverruns,
	pi_l1plistlatestalls
} mypapi_interpretations;

#define MYPAPI_SETS 6
void mypapi_init();
void mypapi_start(int i);
mypapi_counters mypapi_stop();
void mypapi_free();

mypapi_counters mypapi_merge_counters(mypapi_counters *counters1, mypapi_counters *counters2);
void mypapi_print_counters(mypapi_counters *counters);

#define NANO  (1e-9)
#define MICRO (1e-6)
#define MILLI (1e-3)

#define KILO (1e3)
#define MEGA (1e6)
#define GIGA (1e9)
#define TERA (1e12)

#define KIBI (1024.0)
#define MEBI (1024.0*1024.0)
#define GIBI (1024.0*1024.0*1024.0)
#define TEBI (1024.0*1024.0*1024.0*1024.0)

#endif
