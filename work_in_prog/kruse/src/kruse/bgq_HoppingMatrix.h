/*
 * bgq_HoppingMatrix.h
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_HOPPINGMATRIX_H_
#define BGQ_HOPPINGMATRIX_H_

#include "bgq_spinorfield.h"
#include "bgq_utils.h"


#ifndef BGQ_HOPPING_MATRIX_C_
#define EXTERN_INLINE EXTERN_INLINE_DECLARATION
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE EXTERN_INLINE_DEFINITION
#define EXTERN_FIELD
#define EXTERN_INIT(val) = (val)
#endif



void bgq_HoppingMatrix(bool isOdd, bgq_weylfield_controlblock *targetfield, bgq_weylfield_controlblock *spinorfield, bgq_hmflags opts);


typedef struct {
	bool isOdd_src;
	bool isOdd_dst;
	bgq_weylfield_controlblock *targetfield;
	bgq_weylfield_controlblock *spinorfield;
	ucoord ic_begin;
	ucoord ic_end;
	bool noprefetchstream;
} bgq_HoppingMatrix_workload;

void bgq_HoppingMatrix_work(bgq_HoppingMatrix_workload *work,  bool nokamul, bgq_spinorfield_layout layout);


#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_HOPPINGMATRIX_H_ */
