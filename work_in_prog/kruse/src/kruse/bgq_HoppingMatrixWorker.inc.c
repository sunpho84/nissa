
#ifndef BGQ_HOPPINGMATRIXWORKER_INC_
#include "bgq_utils.h"
#include "bgq_HoppingMatrix.h"
#include "bgq_gaugefield.h"

#include "../boundary.h"

#include <stdbool.h>

#define PRECISION double


void bgq_HoppingMatrix_worker(void *arg, size_t tid, size_t threads, bool kamul, bool readFulllayout)
#endif

{
	bgq_HoppingMatrix_workload *work = arg;
	bool isOdd_src = work->isOdd_src;
	bool isOdd_dst = work->isOdd_dst;
	bgq_weylfield_controlblock * restrict spinorfield = work->spinorfield;
	bgq_weylfield_controlblock * restrict targetfield = work->targetfield;
	ucoord ic_begin = work->ic_begin;
	ucoord ic_end = work->ic_end;
	bool noprefetchstream = work->noprefetchstream;

	bgq_vector4double_decl(qka0);
	bgq_complxval_splat(qka0,ka0);
	bgq_vector4double_decl(qka1);
	bgq_complxval_splat(qka1,ka1);
	bgq_vector4double_decl(qka2);
	bgq_complxval_splat(qka2,ka2);
	bgq_vector4double_decl(qka3);
	bgq_complxval_splat(qka3,ka3);

	const size_t workload = ic_end - ic_begin;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = ic_begin + tid*threadload;
	const size_t end = min_sizet(ic_end, begin+threadload);

	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_gaugefield_fromCollapsed[isOdd_src][begin]);
		if (readFulllayout) {
			bgq_prefetch_forward(&spinorfield->BGQ_SEC_FULLLAYOUT[begin]);
		} else {
			bgq_prefetch_forward(&spinorfield->BGQ_SEC_WEYLLAYOUT[begin]);
		}
		bgq_prefetch_forward(&targetfield->BGQ_SENDPTR[isOdd_dst][begin]);
	}

	bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd_src][begin];
	gaugesite = (bgq_gaugesite*)(((uint8_t*)gaugesite)-32);
	bgq_weylsite *weylsite = &spinorfield->BGQ_SEC_WEYLLAYOUT[begin];
	weylsite = (bgq_weylsite*)(((uint8_t*)weylsite) - BGQ_VECTORSIZE);

	for (ucoord ic = begin; ic<end; ic+=1) {
#ifndef NDEBUG
		ucoord ih = bgq_collapsed2halfvolume(isOdd_src, ic);
		ucoord t1 = bgq_halfvolume2t1(isOdd_src, ih);
		ucoord t2 = bgq_halfvolume2t2(isOdd_src, ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);
#endif
		bgq_weyl_ptr_t * restrict targetptrs = &targetfield->BGQ_SENDPTR[isOdd_dst][ic];

		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);

		if (readFulllayout) {
			bgq_su3_matrix_prefetch(gaugesite);
			bgq_spinorsite *spinorsite = &spinorfield->BGQ_SEC_FULLLAYOUT[ic];
			//assert(spinorsite->s[1][0][0]!=0);
			//bgq_su3_spinor_prefetch_double(&spinorfield->sec_fullspinor[ic+1]); // TODO: This prefetch is too early
			//bgq_HoppingMatrix_loadFulllayout(spinor, spinorsite, t1, t2, x, y, z);
			//bgq_su3_spinor_valgen(spinor);
			bgq_su3_spinor_load(spinor, spinorsite);
					bgq_spinorqpx_expect(spinor, t1, t2, x, y, z);

			#define BGQ_COMPUTEWEYL_INC_ 1
			#define BGQ_COMPUTEWEYL_INSERTPREFETCH bgq_su3_spinor_prefetch(weylsite);
			#include "bgq_ComputeWeyl.inc.c"
		} else {
			#define BGQ_READWEYLLAYOUT_INC_ 1
			#define BGQ_READWEYLLAYOUT_INSERTPREFETCH bgq_su3_matrix_prefetch_double(gaugesite);
			#include "bgq_ReadWeyllayout.inc.c"

			#define BGQ_COMPUTEWEYL_INC_ 1
			#define BGQ_COMPUTEWEYL_INSERTPREFETCH bgq_su3_weyl_prefetch(weylsite);
			#include "bgq_ComputeWeyl.inc.c"
		}
	}
}


#undef BGQ_HOPPINGMATRIXWORKER_INC_
