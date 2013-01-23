/*
 * bgq_workers.h
 *
 *  Created on: Nov 19, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_WORKERS_H_
#define BGQ_WORKERS_H_

#include "bgq_spinorfield.h"
#include "bgq_dispatch.h"

#include <string.h>
#include <stdbool.h>


#define BGQ_SPINORFIELD_GENWORKER(name) \
	static void NAME2(name,readFull)(void *arg_untyped, size_t tid, size_t threads) { \
		name(arg_untyped, tid, threads, false, false, false, false); \
	} \
	static void NAME2(name,readWeyl)(void *arg_untyped, size_t tid, size_t threads) { \
		name(arg_untyped, tid, threads, true, false, false, false); \
	} \
	static void NAME2(name,readFullSloppy)(void *arg_untyped, size_t tid, size_t threads) { \
		name(arg_untyped, tid, threads, false, true, false, false); \
	} \
	static void NAME2(name,readWeylSloppy)(void *arg_untyped, size_t tid, size_t threads) { \
		name(arg_untyped, tid, threads, true, true, false, false); \
	} \
	static void NAME2(name,readLegacy)(void *arg_untyped, size_t tid, size_t threads) { \
		name(arg_untyped, tid, threads, false, false, false, true); \
	} \
	\
	bgq_worker_func NAME3(g,name,list)[BGQ_SPINORFIELD_LAYOUT_COUNT] = { \
		&NAME2(name,readFull), &NAME2(name,readWeyl), &NAME2(name,readFullSloppy), &NAME2(name,readWeylSloppy), NULL,NULL,NULL,NULL, NAME2(name,readLegacy) \
	};


#define PRECISION_ISSLOPPY_double false
#define PRECISION_ISSLOPPY_float true
#define PRECISION_ISSLOPPY NAME2(PRECISION_ISSLOPPY,PRECISION)

#define PRECISION_SIZEOF_double 8
#define PRECISION_SIZEOF_float 4
#define PRECISION_SIZEOF NAME2(PRECISION_SIZEOF,PRECISION)

#define PRECISION_COMPLEX_SIZEOF (2*PRECISION_SIZEOF)

typedef struct {
	bool isOdd;
	bgq_weylfield_controlblock *field;
	bgq_hmflags opts;
} bgq_unvectorize_workload;

typedef struct {
	bool isOdd;
	bgq_weylfield_controlblock *field;
	spinor *target;
} bgq_copyToLegacy_workload;

typedef struct {
	bool isOdd;
	spinor *source;
	bgq_weylfield_controlblock *target;
} bgq_copyFromLegacy_workload;

typedef struct {
	bgq_weylfield_controlblock *field;
	bool isOdd;
} bgq_spinorfield_rewrite_work;

void bgq_HoppingMatrix_unvectorize_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_HoppingMatrix_unvectorize_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_HoppingMatrix_unvectorize NAME2(bgq_HoppingMatrix_unvectorize,PRECISION)

void bgq_HoppingMatrix_worker_datamove_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_HoppingMatrix_worker_datamove_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_HoppingMatrix_worker_datamove NAME2(bgq_HoppingMatrix_worker_datamove,PRECISION)

void bgq_HoppingMatrix_datamovet_worker_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_HoppingMatrix_datamovet_worker_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_HoppingMatrix_datamovet_worker NAME2(bgq_HoppingMatrix_datamovet_worker,PRECISION)


#if 0
void bgq_copyToLegacy_worker_fulllayout_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_copyToLegacy_worker_fulllayout_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_copyToLegacy_worker_fulllayout NAME2(bgq_copyToLegacy_worker_fulllayout,PRECISION)

void bgq_copyToLegacy_worker_weyllayout_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_copyToLegacy_worker_weyllayout_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_copyToLegacy_worker_weyllayout NAME2(bgq_copyToLegacy_worker_weyllayout,PRECISION)
#endif


void bgq_copyFromLegacy_worker_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_copyFromLegacy_worker_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_copyFromLegacy_worker NAME2(bgq_copyFromLegacy_worker,PRECISION)

extern bgq_worker_func g_bgq_spinorfield_rewrite_worker_double_list[BGQ_SPINORFIELD_LAYOUT_COUNT];
extern bgq_worker_func g_bgq_spinorfield_rewrite_worker_float_list[BGQ_SPINORFIELD_LAYOUT_COUNT];
#define bgq_spinorfield_rewrite_worker NAME2(bgq_spinorfield_rewrite_worker,PRECISION)

#endif /* BGQ_WORKERS_H_ */
