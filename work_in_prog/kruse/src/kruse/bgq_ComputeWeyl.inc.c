/*
 * IBM xlc refuses to inline functions that big, therefore we have to include the function using the preprocessor
 */


#ifndef BGQ_COMPUTEWEYL_INC_
#include "bgq_qpx.h"
#include "bgq_spinorfield.h"
#include "bgq_gaugefield.h"

#include <stdbool.h>

#define PRECISION double

void bgq_HoppingMatrix_compute_storeWeyllayout_raw(bgq_weyl_ptr_t *targetptrs, bgq_gaugesite *gaugesite, bgq_su3_spinor_params(spinor), bgq_params(qka0), bgq_params(qka1), bgq_params(qka2), bgq_params(qka3), ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z, bool kamul)
#endif

#ifndef BGQ_COMPUTEWEYL_INSERTPREFETCH
#define BGQ_COMPUTEWEYL_INSERTPREFETCH
#endif

{
#ifndef NDEBUG
	bool isOdd_src = bgq_local2isOdd(t1,x,y,z);
	assert(isOdd_src == bgq_local2isOdd(t2,x,y,z));
	bool isOdd_dst = !isOdd_src;

	ucoord ic_src = bgq_local2collapsed(t1, x, y, z);
	assert(ic_src == bgq_local2collapsed(t1, x, y, z));
#endif

	REORDER_BARRIER // These keep the compiler from spilling because of excessive instruction reordering
	bgq_prefetch(&targetptrs->d[TUP]);
	bgq_su3_matrixnext_prefetch(gaugesite);
	//bgq_su3_weyl_prestore_double(targetptrs->d[TUP]);

	// T+ /////////////////////////////////////////////////////////////////////////
	{
				bgq_setdesc(BGQREF_TUP_SOURCE,"BGQREF_TUP_SOURCE");
				bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP_SOURCE, bgq_cmplxval1(spinor_v0_c0));
				bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP_SOURCE, bgq_cmplxval2(spinor_v0_c0));

		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_reduce_weyl_tup(weyl_tdown, spinor);
				bgq_setdesc(BGQREF_TUP, "BGQREF_TUP");
				bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP, bgq_cmplxval1(weyl_tdown_v0_c0));
				bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP, bgq_cmplxval2(weyl_tdown_v0_c0));

		bgq_su3_mdecl(gauge_tup);
		bgq_qvlfduxa(gauge_tup_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c22, gaugesite, 32);
				bgq_gaugeqpx_expect(gauge_tup, t1, t2, x, y, z, TUP, true);
				bgq_setdesc(BGQREF_TUP_GAUGE, "BGQREF_TUP_GAUGE");
				bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP_GAUGE, bgq_cmplxval1(gauge_tup_c00));
				bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP_GAUGE, bgq_cmplxval2(gauge_tup_c00));

		bgq_su3_weyl_mvmul(weyl_tdown, gauge_tup, weyl_tdown);
				bgq_setdesc(BGQREF_TUP_WEYL,"BGQREF_TUP_WEYL");
				bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP_WEYL, bgq_cmplxval1(weyl_tdown_v0_c0));
				bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP_WEYL, bgq_cmplxval2(weyl_tdown_v0_c0));
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_tdown, qka0, weyl_tdown);
		}
				bgq_setdesc(BGQREF_TUP_KAMUL,"BGQREF_TUP_KAMUL");
				bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP_KAMUL, bgq_cmplxval1(weyl_tdown_v0_c0));
				bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP_KAMUL, bgq_cmplxval2(weyl_tdown_v0_c0));

		//bgq_su3_weyl_zeroload(targetptrs->d[TUP]);
		//bgq_su3_weyl_prestore_double(targetptrs->d[TUP]);
		bgq_su3_weyl_store(targetptrs->d[TUP], weyl_tdown);
				bgq_weylvec_written(targetptrs->d[TUP], t1, t2, x, y, z, TDOWN, true);
	}

	REORDER_BARRIER
	bgq_prefetch(&targetptrs->d[TDOWN]); // Note: just one of these should be necessary because the complete bgq_weyl_ptr_t fits into one 64-byte cache line -- but performance is better with them; maybe they are thrown out of cache prematurely
	bgq_su3_matrixnext_prefetch(gaugesite);
	//bgq_su3_weyl_prestore_double(targetptrs->d[TDOWN]);

	// T- /////////////////////////////////////////////////////////////////////////
	{
				bgq_setdesc(BGQREF_TDOWN_SOURCE, "BGQREF_TDOWN_SOURCE");
				bgq_setbgqvalue_src(t1, x, y, z, TUP, BGQREF_TDOWN_SOURCE, bgq_cmplxval1(spinor_v0_c0));
				bgq_setbgqvalue_src(t2, x, y, z, TUP, BGQREF_TDOWN_SOURCE, bgq_cmplxval2(spinor_v0_c0));

		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_reduce_weyl_tdown(weyl_tup, spinor);
				bgq_setdesc(BGQREF_TDOWN, "BGQREF_TDOWN");
				bgq_setbgqvalue_src(t1, x, y, z, TUP, BGQREF_TDOWN, bgq_cmplxval1(weyl_tup_v0_c2));
				bgq_setbgqvalue_src(t2, x, y, z, TUP, BGQREF_TDOWN, bgq_cmplxval2(weyl_tup_v0_c2));

		bgq_su3_mdecl(gauge_tdown);
		bgq_qvlfduxa(gauge_tdown_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c22, gaugesite, 32);
				bgq_gaugeqpx_expect(gauge_tdown, t1, t2, x, y, z, TDOWN, true);
				bgq_setdesc(BGQREF_TDOWN_GAUGE, "BGQREF_TDOWN_GAUGE");
				bgq_setbgqvalue_src(t1, x, y, z, TUP, BGQREF_TDOWN_GAUGE, bgq_cmplxval1(gauge_tdown_c00));
				bgq_setbgqvalue_src(t2, x, y, z, TUP, BGQREF_TDOWN_GAUGE, bgq_cmplxval2(gauge_tdown_c00));

		bgq_su3_weyl_mvinvmul(weyl_tup, gauge_tdown, weyl_tup);
				bgq_setdesc(BGQREF_TDOWN_WEYL,"BGQREF_TDOWN_WEYL");
				bgq_setbgqvalue_src(t1, x, y, z, TUP, BGQREF_TDOWN_WEYL, bgq_cmplxval1(weyl_tup_v0_c0));
				bgq_setbgqvalue_src(t2, x, y, z, TUP, BGQREF_TDOWN_WEYL, bgq_cmplxval2(weyl_tup_v0_c0));
		if (kamul) {
			bgq_su3_weyl_cjgmul(weyl_tup, qka0, weyl_tup);
		}
				bgq_setdesc(BGQREF_TDOWN_KAMUL,"BGQREF_TDOWN_KAMUL");
				bgq_setbgqvalue_src(t1, x, y, z, TUP, BGQREF_TDOWN_KAMUL, bgq_cmplxval1(weyl_tup_v0_c0));
				bgq_setbgqvalue_src(t2, x, y, z, TUP, BGQREF_TDOWN_KAMUL, bgq_cmplxval2(weyl_tup_v0_c0));

		//bgq_su3_weyl_zeroload(targetptrs->d[TDOWN]); //TODO: remove?
		//bgq_su3_weyl_prestore_double(targetptrs->d[TDOWN]);
		bgq_su3_weyl_store(targetptrs->d[TDOWN], weyl_tup);
				bgq_weylvec_written(targetptrs->d[TDOWN], t1, t2, x, y, z, TUP, true);
	}

	REORDER_BARRIER
	bgq_prefetch(&targetptrs->d[XUP]);
	bgq_su3_matrixnext_prefetch(gaugesite);
	//bgq_su3_weyl_prestore_double(targetptrs->d[XUP]);

	// X+ /////////////////////////////////////////////////////////////////////////
	{
			bgq_su3_weyl_decl(weyl_xdown);
			bgq_su3_reduce_weyl_xup(weyl_xdown, spinor);
					bgq_setdesc(BGQREF_XUP, "BGQREF_XUP");
					bgq_setbgqvalue_src(t1, x, y, z, XDOWN, BGQREF_XUP, bgq_cmplxval1(weyl_xdown_v0_c0));
					bgq_setbgqvalue_src(t2, x, y, z, XDOWN, BGQREF_XUP, bgq_cmplxval2(weyl_xdown_v0_c0));

			bgq_su3_mdecl(gauge_xup);
			bgq_qvlfduxa(gauge_xup_c00, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c01, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c02, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c10, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c11, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c12, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c20, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c21, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c22, gaugesite, 32);
					bgq_gaugeqpx_expect(gauge_xup, t1, t2, x, y, z, XUP, true);
					bgq_setdesc(BGQREF_XUP_GAUGE, "BGQREF_XUP_GAUGE");
					bgq_setbgqvalue_src(t1, x, y, z, XDOWN, BGQREF_XUP_GAUGE, bgq_cmplxval1(gauge_xup_c00));
					bgq_setbgqvalue_src(t2, x, y, z, XDOWN, BGQREF_XUP_GAUGE, bgq_cmplxval2(gauge_xup_c00));
			bgq_su3_weyl_mvmul(weyl_xdown, gauge_xup, weyl_xdown);
					bgq_setdesc(BGQREF_XUP_WEYL,"BGQREF_XUP_WEYL");
					bgq_setbgqvalue_src(t1, x, y, z, XDOWN, BGQREF_XUP_WEYL, bgq_cmplxval1(weyl_xdown_v0_c0));
					bgq_setbgqvalue_src(t2, x, y, z, XDOWN, BGQREF_XUP_WEYL, bgq_cmplxval2(weyl_xdown_v0_c0));
			if (kamul) {
				bgq_su3_weyl_cmul(weyl_xdown, qka1, weyl_xdown);
			}
					bgq_setdesc(BGQREF_XUP_KAMUL,"BGQREF_XUP_KAMUL");
					bgq_setbgqvalue_src(t1, x, y, z, XDOWN, BGQREF_XUP_KAMUL, bgq_cmplxval1(weyl_xdown_v0_c0));
					bgq_setbgqvalue_src(t2, x, y, z, XDOWN, BGQREF_XUP_KAMUL, bgq_cmplxval2(weyl_xdown_v0_c0));

			//bgq_su3_weyl_zeroload(targetptrs->d[XUP]);
			//bgq_su3_weyl_prestore_double(targetptrs->d[XUP]);
			bgq_su3_weyl_store(targetptrs->d[XUP], weyl_xdown);
					bgq_weylvec_written(targetptrs->d[XUP], t1, t2, x,y,z,XDOWN, true);
		}

	REORDER_BARRIER
	bgq_prefetch(&targetptrs->d[XDOWN]);
	bgq_su3_matrixnext_prefetch(gaugesite);
	//bgq_su3_weyl_prestore_double(targetptrs->d[XDOWN]);

	// X- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_reduce_weyl_xdown(weyl_xup, spinor);
				bgq_setdesc(BGQREF_XDOWN, "BGQREF_XDOWN");
				bgq_setbgqvalue_src(t1, x, y, z, XUP, BGQREF_XDOWN, bgq_cmplxval1(weyl_xup_v0_c0));
				bgq_setbgqvalue_src(t2, x, y, z, XUP, BGQREF_XDOWN, bgq_cmplxval2(weyl_xup_v0_c0));

		bgq_su3_mdecl(gauge_xdown);
		bgq_qvlfduxa(gauge_xdown_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c22, gaugesite, 32);
				bgq_gaugeqpx_expect(gauge_xdown, t1, t2, x, y, z, XDOWN, true);
				bgq_setdesc(BGQREF_XDOWN_GAUGE, "BGQREF_XDOWN_GAUGE");
				bgq_setbgqvalue_src(t1, x, y, z, XUP, BGQREF_XDOWN_GAUGE, bgq_cmplxval1(gauge_xdown_c00));
				bgq_setbgqvalue_src(t2, x, y, z, XUP, BGQREF_XDOWN_GAUGE, bgq_cmplxval2(gauge_xdown_c00));

		bgq_su3_weyl_mvinvmul(weyl_xup, gauge_xdown, weyl_xup);
				bgq_setdesc(BGQREF_XDOWN_WEYL, "BGQREF_XDOWN_WEYL");
				bgq_setbgqvalue_src(t1, x, y, z, XUP, BGQREF_XDOWN_WEYL, bgq_cmplxval1(weyl_xup_v0_c0));
				bgq_setbgqvalue_src(t2, x, y, z, XUP, BGQREF_XDOWN_WEYL, bgq_cmplxval2(weyl_xup_v0_c0));
		if (kamul) {
			bgq_su3_weyl_cjgmul(weyl_xup, qka1, weyl_xup);
		}
				bgq_setdesc(BGQREF_XDOWN_KAMUL,"BGQREF_XDOWN_KAMUL");
				bgq_setbgqvalue_src(t1, x, y, z, XUP, BGQREF_XDOWN_KAMUL, bgq_cmplxval1(weyl_xup_v0_c0));
				bgq_setbgqvalue_src(t2, x, y, z, XUP, BGQREF_XDOWN_KAMUL, bgq_cmplxval2(weyl_xup_v0_c0));

		//bgq_su3_weyl_zeroload(targetptrs->d[XDOWN]);
		//bgq_su3_weyl_prestore_double(targetptrs->d[ZDOWN]);
		bgq_su3_weyl_store(targetptrs->d[XDOWN], weyl_xup);
				bgq_weylvec_written(targetptrs->d[XDOWN], t1, t2, x,y,z,XUP, true);
	}

	REORDER_BARRIER
	bgq_prefetch(&targetptrs->d[YUP]);
	bgq_su3_matrixnext_prefetch(gaugesite);
	//bgq_su3_weyl_prestore_double(targetptrs->d[YUP]);

	// Y+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_reduce_weyl_yup(weyl_ydown, spinor);

		bgq_su3_mdecl(gauge_yup);
		bgq_qvlfduxa(gauge_yup_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c22, gaugesite, 32);
				bgq_gaugeqpx_expect(gauge_yup, t1, t2, x, y, z, YUP, true);

		bgq_su3_weyl_mvmul(weyl_ydown, gauge_yup, weyl_ydown);
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_ydown, qka2, weyl_ydown);
		}

		//bgq_su3_weyl_zeroload(targetptrs->d[YUP]);
		//bgq_su3_weyl_prestore_double(targetptrs->d[YUP]);
		bgq_su3_weyl_store(targetptrs->d[YUP], weyl_ydown);
				bgq_weylvec_written(targetptrs->d[YUP], t1, t2, x,y,z,YDOWN, true);
	}

	REORDER_BARRIER
	bgq_prefetch(&targetptrs->d[YDOWN]);
	bgq_su3_matrixnext_prefetch(gaugesite);
	//bgq_su3_weyl_prestore_double(targetptrs->d[YDOWN]);

	// Y- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_reduce_weyl_ydown(weyl_yup, spinor);

		bgq_su3_mdecl(gauge_ydown);
		bgq_qvlfduxa(gauge_ydown_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c22, gaugesite, 32);
				bgq_gaugeqpx_expect(gauge_ydown, t1, t2, x, y, z, YDOWN, true);

		bgq_su3_weyl_mvinvmul(weyl_yup, gauge_ydown, weyl_yup);
		if (kamul) {
			bgq_su3_weyl_cjgmul(weyl_yup, qka2, weyl_yup);
		}

		//bgq_su3_weyl_zeroload(targetptrs->d[YDOWN]);
		//bgq_su3_weyl_prestore_double(targetptrs->d[YDOWN]);
		bgq_su3_weyl_store(targetptrs->d[YDOWN], weyl_yup);
				bgq_weylvec_written(targetptrs->d[YDOWN], t1, t2, x,y,z,YUP, true);
	}

	REORDER_BARRIER
	bgq_prefetch(&targetptrs->d[ZUP]);
	bgq_su3_matrixnext_prefetch(gaugesite);
	//bgq_su3_weyl_prestore_double(targetptrs->d[ZUP]);

	// Z+ /////////////////////////////////////////////////////////////////////////
	{
			bgq_su3_weyl_decl(weyl_zdown);
			bgq_su3_reduce_weyl_zup(weyl_zdown, spinor);

			bgq_su3_mdecl(gauge_zup);
			bgq_qvlfduxa(gauge_zup_c00, gaugesite, 32);
			bgq_qvlfduxa(gauge_zup_c01, gaugesite, 32);
			bgq_qvlfduxa(gauge_zup_c02, gaugesite, 32);
			bgq_qvlfduxa(gauge_zup_c10, gaugesite, 32);
			bgq_qvlfduxa(gauge_zup_c11, gaugesite, 32);
			bgq_qvlfduxa(gauge_zup_c12, gaugesite, 32);
			bgq_qvlfduxa(gauge_zup_c20, gaugesite, 32);
			bgq_qvlfduxa(gauge_zup_c21, gaugesite, 32);
			bgq_qvlfduxa(gauge_zup_c22, gaugesite, 32);
					bgq_gaugeqpx_expect(gauge_zup, t1, t2, x, y, z, ZUP, true);

			bgq_su3_weyl_mvmul(weyl_zdown, gauge_zup, weyl_zdown);
			if (kamul) {
				bgq_su3_weyl_cmul(weyl_zdown, qka3, weyl_zdown);
			}

			//bgq_su3_weyl_zeroload(targetptrs->d[ZUP]);
			//bgq_su3_weyl_prestore_double(targetptrs->d[ZUP]);
			bgq_su3_weyl_store(targetptrs->d[ZUP], weyl_zdown);
					bgq_weylvec_written(targetptrs->d[ZUP], t1, t2, x,y,z, ZDOWN, true);
		}

	REORDER_BARRIER
	bgq_prefetch(&targetptrs->d[ZDOWN]);
	//bgq_su3_weyl_prestore_double(targetptrs->d[ZDOWN]);
	BGQ_COMPUTEWEYL_INSERTPREFETCH

	// Z- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_reduce_weyl_zdown(weyl_zup, spinor);

		bgq_su3_mdecl(gauge_zdown);
		bgq_qvlfduxa(gauge_zdown_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c22, gaugesite, 32);
				bgq_gaugeqpx_expect(gauge_zdown, t1, t2, x, y, z, ZDOWN, true);

		bgq_su3_weyl_mvinvmul(weyl_zup, gauge_zdown, weyl_zup);
		if (kamul) {
			bgq_su3_weyl_cjgmul(weyl_zup, qka3, weyl_zup);
		}

		//bgq_su3_weyl_zeroload(targetptrs->d[ZDOWN]);
		//bgq_su3_weyl_prestore_double(targetptrs->d[ZDOWN]);
		bgq_su3_weyl_store(targetptrs->d[ZDOWN], weyl_zup);
				bgq_weylvec_written(targetptrs->d[ZDOWN], t1, t2, x,y,z,ZUP, true);
	}

	REORDER_BARRIER
}

#undef BGQ_COMPUTEWEYL_INSERTPREFETCH
#undef BGQ_COMPUTEWEYL_INC_
