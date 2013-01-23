/*
 * bgq_HoppingMatrix.c
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */

#define BGQ_HOPPING_MATRIX_C_
#include "bgq_HoppingMatrix.h"

#include "bgq_field.h"
#include "bgq_spinorfield.h"
#include "bgq_qpx.h"
#include "bgq_dispatch.h"
#include "bgq_comm.h"
#include "bgq_workers.h"
#include "bgq_gaugefield.h"

#include "../interface/interface.h"
//#include "../boundary.h"
//#include "../update_backward_gauge.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>


static void bgq_HoppingMatrix_nokamul_worker_readFulllayout_double(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,false,true);

	bool const kamul = false;
	bool const readFulllayout = true;

#define PRECISION double
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_nokamul_worker_readFulllayout_float(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,false,true);

	bool const kamul = false;
	bool const readFulllayout = true;

#define PRECISION float
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_kamul_worker_readFulllayout_double(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,true,true);

	bool const kamul = true;
	bool const readFulllayout = true;

#define PRECISION double
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_kamul_worker_readFulllayout_float(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,true,true);

	bool const kamul = true;
	bool const readFulllayout = true;

#define PRECISION float
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_nokamul_worker_readWeyllayout_double(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,false,false);

	bool const kamul = false;
	bool const readFulllayout = false;

#define PRECISION double
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_nokamul_worker_readWeyllayout_float(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,false,false);

	bool const kamul = false;
	bool const readFulllayout = false;

#define PRECISION float
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_kamul_worker_readWeyllayout_double(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,true,false);

	bool const kamul = true;
	bool const readFulllayout = false;

#define PRECISION double
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_kamul_worker_readWeyllayout_float(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,true,false);

	bool const kamul = true;
	bool const readFulllayout = false;

#define PRECISION float
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}


void bgq_HoppingMatrix_work(bgq_HoppingMatrix_workload *work,  bool nokamul, bgq_spinorfield_layout layout) {
	assert(work);
	size_t sites = work->ic_end - work->ic_begin;
	uint64_t old = flopaccumulator;

	bgq_worker_func func;
	if (nokamul) {
		switch (layout) {
		case ly_full_double:
			func = &bgq_HoppingMatrix_nokamul_worker_readFulllayout_double;
			break;
		case ly_full_float:
			func = &bgq_HoppingMatrix_nokamul_worker_readFulllayout_float;
			break;
		case ly_weyl_double:
			func = &bgq_HoppingMatrix_nokamul_worker_readWeyllayout_double;
			break;
		case ly_weyl_float:
			func = &bgq_HoppingMatrix_nokamul_worker_readWeyllayout_float;
			break;
		default:
			UNREACHABLE
		}
	} else {
		switch (layout) {
		case ly_full_double:
			func = &bgq_HoppingMatrix_kamul_worker_readFulllayout_double;
			break;
		case ly_full_float:
			func = &bgq_HoppingMatrix_kamul_worker_readFulllayout_float;
			break;
		case ly_weyl_double:
			func = &bgq_HoppingMatrix_kamul_worker_readWeyllayout_double;
			break;
		case ly_weyl_float:
			func = &bgq_HoppingMatrix_kamul_worker_readWeyllayout_float;
			break;
		default:
			UNREACHABLE
		}
	}
	bgq_master_call(func, work);
	uint64_t flopPerSite = (
			/*accum spinor*/	7/*dirs*/ * (4 * 3)/*cmplx per spinor*/ * 2/*flops accum*/ +
			/*su3 mul*/			8/*dirs*/ * 2/*su3vec per weyl*/ * (9*6/*cmpl mul*/ + 6*2/*cmplx add*/)/*flop per su3 mv-mul*/
	);
	if (!nokamul)
		flopPerSite += /*kamul*/			8/*dirs*/ * (2 * 3)/*cmplx per weyl*/ * 6/*flops cmplx mul*/;
	if (layout & ly_weyl)
		flopPerSite += /*weyl reduce*/      8/*dirs*/ * (2 * 3)/*cmplx per weyl*/ * 2/*flops*/;
	flopaccumulator += sites * PHYSICAL_LK * flopPerSite;
}


void bgq_HoppingMatrix(bool isOdd, bgq_weylfield_controlblock *targetfield, bgq_weylfield_controlblock *spinorfield, bgq_hmflags opts) {
	assert(spinorfield);
	assert(targetfield);
	assert(spinorfield != targetfield);

	bool nocomm = opts & hm_nocom;
	bool nooverlap = opts & hm_nooverlap;
	bool nokamul = opts & hm_nokamul;
	bool nodistribute = opts & hm_nodistribute;
	bool nodatamove = opts & hm_nodatamove;
	bool nobody = opts & hm_nobody;
	bool nospi = opts & hm_nospi;
	bool noprefetchstream = opts & hm_noprefetchstream;
	bool floatprecision = opts & hm_floatprecision;

	bgq_spinorfield_layout layout = bgq_spinorfield_prepareRead(spinorfield, !isOdd, true, !floatprecision, floatprecision, false, false);
	bgq_spinorfield_prepareWrite(targetfield, isOdd, floatprecision ? ly_weyl_float : ly_weyl_double, targetfield==spinorfield);


	if(g_update_gauge_copy) {
		update_backward_gauge(g_gauge_field);
	}


	// 0. Expect data from other neighbor node
	if (!nocomm) {
		bgq_comm_recv(nospi, floatprecision, targetfield);
	}


	// 1. Distribute
	// Compute surface and put data into the send buffers
	if ((PHYSICAL_SURFACE > 0) && !nodistribute) {
		bgq_master_sync();
		static bgq_HoppingMatrix_workload work_surface;
		work_surface.isOdd_src = !isOdd;
		work_surface.isOdd_dst = isOdd;
		work_surface.targetfield = targetfield;
		work_surface.spinorfield = spinorfield;
		work_surface.ic_begin = bgq_surface2collapsed(0);
		work_surface.ic_end = bgq_surface2collapsed(PHYSICAL_SURFACE-1)+1;
		work_surface.noprefetchstream = noprefetchstream;
		bgq_HoppingMatrix_work(&work_surface, nokamul, layout);
	}


	if (BGQ_UNVECTORIZE && COMM_T && !nodatamove) {
		bgq_master_sync();
		static bgq_unvectorize_workload work_unvectorize;
		work_unvectorize.isOdd = isOdd;
		work_unvectorize.field = targetfield;
		work_unvectorize.opts = opts;
		if (floatprecision)
			bgq_master_call(&bgq_HoppingMatrix_unvectorize_float, &work_unvectorize);
		else
			bgq_master_call(&bgq_HoppingMatrix_unvectorize_double, &work_unvectorize);
	}


// 2. Start communication
	if ((PHYSICAL_SURFACE > 0) && !nocomm) {
		bgq_master_sync(); // Wait for threads to finish surface before sending it
		//TODO: ensure there are no other communications pending
		bgq_comm_send();
		//targetfield->waitingForRecv = true;
		if (!nodatamove) {
			targetfield->pendingDatamove = true;
		}
	}

	targetfield->hmflags = opts;

	if (nooverlap) {
		// Do not wait until data is required, but do it here
		bgq_spinorfield_prepareRead(targetfield, isOdd, true, true, true, true, false);
		//bgq_spinorfield_setup(targetfield, isOdd, false, false, true, false, false);
	}


// 3. Compute the body
	if ((PHYSICAL_BODY > 0) && !nobody) {
		bgq_master_sync();
		static bgq_HoppingMatrix_workload work_body;
		work_body.isOdd_src = !isOdd;
		work_body.isOdd_dst = isOdd;
		work_body.targetfield = targetfield;
		work_body.spinorfield = spinorfield;
		work_body.ic_begin = bgq_body2collapsed(0);
		work_body.ic_end = bgq_body2collapsed(PHYSICAL_BODY-1)+1;
		work_body.noprefetchstream = noprefetchstream;
		bgq_HoppingMatrix_work(&work_body, nokamul, layout);
		if (!COMM_T && !nodatamove) {
			// Copy the data from HALO_T into the required locations
			bgq_master_sync();
			static bgq_work_datamove work_datamovet;
			work_datamovet.spinorfield = targetfield;
			work_datamovet.opts = opts;
			if (floatprecision)
				bgq_master_call(&bgq_HoppingMatrix_datamovet_worker_float, &work_datamovet);
			else
				bgq_master_call(&bgq_HoppingMatrix_datamovet_worker_double, &work_datamovet);
		}
	}


// 4. Wait for the communication to finish
	/* Defer to bgq_spinorfield_setup as soon as the data is actually required */


// 5. Move received to correct location
	/* Done in bgq_spinorfield_setup whoever is using the field next*/


// 6. Compute the surface
	/* Done by funcs calling readWeyllayout */
}

