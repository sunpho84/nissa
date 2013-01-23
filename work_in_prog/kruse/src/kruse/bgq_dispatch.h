/*
 * bgq_dispatch.h
 *
 *  Created on: Oct 15, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_DISPATCH_H_
#define BGQ_DISPATCH_H_

#include "bgq_utils.h"

#include <string.h>
#include <stdbool.h>

#ifndef BGQ_DISPATCH_C_
#define EXTERN_INLINE EXTERN_INLINE_DECLARATION
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE EXTERN_INLINE_DEFINITION
#define EXTERN_FIELD
#define EXTERN_INIT(val) = (val)
#endif


typedef void (*bgq_worker_func)(void *arg, size_t tid, size_t threads);
typedef int (*bgq_master_func)(void *arg);
typedef int (*bgq_mainlike_func)(int argc, char *argv[]);



EXTERN_FIELD bool g_bgq_dispatch_inparallel EXTERN_INIT(false);
EXTERN_FIELD int g_bgq_dispatch_threads;
int bgq_parallel(bgq_master_func master_func, void *master_arg);
int bgq_parallel_mainlike(bgq_mainlike_func func, int argc, char *argv[]);

void bgq_worker();
void bgq_master_call(bgq_worker_func worker_func, void *arg);
void bgq_master_sync();

void bgq_master_memzero(void *ptr, size_t size);
void bgq_master_memcpy(void *ptrDst, void *ptrSrc, size_t size);

#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_DISPATCH_H_ */
