/*
 * bgq_workers.inc.c
 *
 *  Created on: Nov 19, 2012
 *      Author: meinersbur
 */

#include "bgq_workers.h"

#include "bgq_spinorfield.h"
#include "bgq_field.h"
#include "bgq_comm.h"
#include "bgq_qpx.h"
#include "bgq_utils.h"
#include "bgq_dispatch.h"

#ifndef BGQ_WORKERS_INC_
#define PRECISION double
#undef bgq_HoppingMatrix_unvectorize
#undef bgq_HoppingMatrix_worker_datamove
#undef bgq_HoppingMatrix_datamovet_worker
#undef bgq_copyFromLegacy_worker
#undef bgq_spinorfield_rewrite_worker
#endif



void bgq_HoppingMatrix_unvectorize(void *arg_untyped, size_t tid, size_t threads) {
	assert(BGQ_UNVECTORIZE);
	assert(COMM_T);
	bgq_unvectorize_workload *arg = arg_untyped;
	bool isOdd = arg->isOdd;
	bgq_weylfield_controlblock *field = arg->field;
	bgq_hmflags opts = arg->opts;

	bool noprefetchstream = opts & hm_noprefetchstream;

	const size_t workload_tdown = LOCAL_HALO_T/(PHYSICAL_LP*PHYSICAL_LK);
	const size_t workload_tup = workload_tdown;
	const size_t workload = workload_tdown + workload_tup;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);

	for (size_t i = begin; i < end; ) {
		WORKLOAD_DECL(i, workload);

		if (WORKLOAD_SPLIT(workload_tup)) {
			const size_t beginj = WORKLOAD_PARAM(workload_tup);
			const size_t endj = min_sizet(workload_tup, beginj+threadload);
			if (!noprefetchstream) {
				bgq_prefetch_forward(&g_bgq_sec_temp_tup[2*beginj]);
			}
			for (size_t j = beginj; j < endj; j+=1) {
#ifndef NDEBUG
				size_t offset1 = bgq_pointer2offset(field, &g_bgq_sec_temp_tup[2*j]);
				ucoord index1 = bgq_offset2index(offset1);
				ucoord ic1 = bgq_index2collapsed(isOdd, index1, 0);
				ucoord ih1 = bgq_collapsed2halfvolume(isOdd, ic1);
				ucoord t1 = bgq_halfvolume2t(isOdd, ih1, 0);
				ucoord x1 = bgq_halfvolume2x(ih1);
				ucoord y1 = bgq_halfvolume2y(ih1);
				ucoord z1 = bgq_halfvolume2z(ih1);
				size_t offset2 = bgq_pointer2offset(field, &g_bgq_sec_temp_tup[2*j+1]);
				ucoord index2 = bgq_offset2index(offset2);
				ucoord ic2 = bgq_index2collapsed(isOdd, index2, 0);
				ucoord ih2 = bgq_collapsed2halfvolume(isOdd, ic2);
				ucoord t2 = bgq_halfvolume2t(isOdd, ih2, 0);
				ucoord x2 = bgq_halfvolume2x(ih2);
				ucoord y2 = bgq_halfvolume2y(ih2);
				ucoord z2 = bgq_halfvolume2z(ih2);
#endif

				bgq_su3_weyl_decl(weyl1);
				bgq_su3_weylnextnext_prefetch(&g_bgq_sec_temp_tup[2*j]);
				bgq_su3_weyl_load(weyl1, &g_bgq_sec_temp_tup[2*j]);
						bgq_weylqpxk_expect(weyl1, 1, t1, x1, y1, z1, TDOWN, false);

				bgq_su3_weyl_decl(weyl2);
				bgq_su3_weylnextnext_prefetch(&g_bgq_sec_temp_tup[2*j+1]);
				bgq_su3_weyl_load(weyl2, &g_bgq_sec_temp_tup[2*j+1]);
						bgq_weylqpxk_expect(weyl2, 1, t2, x2, y2, z2, TDOWN, false);

				bgq_su3_weyl_decl(weyl);
				bgq_su3_weyl_rmerge(weyl, weyl1, weyl2);
						bgq_weylqpxk_expect(weyl, 0, t1, x1, y1, z1, TDOWN, false);
						bgq_weylqpxk_expect(weyl, 1, t2, x2, y2, z2, TDOWN, false);

				bgq_su3_weyl_store(&g_bgq_sec_send[TUP][j], weyl);
			}
			i += (endj - beginj);
		} else if (WORKLOAD_SPLIT(workload_tdown)) {
			const size_t beginj = WORKLOAD_PARAM(workload_tdown);
			const size_t endj = min_sizet(workload_tup, beginj+threadload);
			if (!noprefetchstream) {
				bgq_prefetch_forward(&g_bgq_sec_temp_tdown[2*beginj]);
			}
			for (size_t j = beginj; j < endj; j+=1) {
#ifndef NDEBUG
				size_t offset1 = bgq_pointer2offset(field, &g_bgq_sec_temp_tdown[2*j]);
				ucoord index1 = bgq_offset2index(offset1);
				ucoord ic1 = bgq_index2collapsed(isOdd, index1, 0);
				ucoord ih1 = bgq_collapsed2halfvolume(isOdd, ic1);
				ucoord t1 = bgq_halfvolume2t(isOdd, ih1, 1);
				ucoord x1 = bgq_halfvolume2x(ih1);
				ucoord y1 = bgq_halfvolume2y(ih1);
				ucoord z1 = bgq_halfvolume2z(ih1);
				size_t offset2 = bgq_pointer2offset(field, &g_bgq_sec_temp_tdown[2*j+1]);
				ucoord index2 = bgq_offset2index(offset2);
				ucoord ic2 = bgq_index2collapsed(isOdd, index2, 0);
				ucoord ih2 = bgq_collapsed2halfvolume(isOdd, ic2);
				ucoord t2 = bgq_halfvolume2t(isOdd, ih2, 1);
				ucoord x2 = bgq_halfvolume2x(ih2);
				ucoord y2 = bgq_halfvolume2y(ih2);
				ucoord z2 = bgq_halfvolume2z(ih2);
#endif

				bgq_su3_weyl_decl(weyl1);
				bgq_su3_weylnextnext_prefetch(&g_bgq_sec_temp_tdown[2*j]);
				bgq_su3_weyl_load(weyl1, &g_bgq_sec_temp_tdown[2*j]);
						bgq_weylqpxk_expect(weyl1, 0, t1, x1, y1, z1, TUP, false);

				bgq_su3_weyl_decl(weyl2);
				bgq_su3_weylnextnext_prefetch(&g_bgq_sec_temp_tdown[2*j+1]);
				bgq_su3_weyl_load(weyl2, &g_bgq_sec_temp_tdown[2*j+1]);
						bgq_weylqpxk_expect(weyl2, 0, t2, x2, y2, z2, TUP, false);

				bgq_su3_weyl_decl(weyl);
				bgq_su3_weyl_lmerge(weyl, weyl1, weyl2);
						bgq_weylqpxk_expect(weyl, 0, t1, x1, y1, z1, TUP, false);
						bgq_weylqpxk_expect(weyl, 1, t2, x2, y2, z2, TUP, false);

				bgq_su3_weyl_store(&g_bgq_sec_send[TDOWN][j], weyl);
			}
			i += (endj - beginj);
		} else {
			UNREACHABLE
		}

		WORKLOAD_CHECK
	}
}


static inline void bgq_HoppingMatrix_worker_datamove_recvxyz(bgq_weylfield_controlblock *spinorfield, bgq_direction d, bool isOdd, size_t beginj, size_t endj, bool noprefetchstream) {
	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_sec_recv[d][beginj]);
		bgq_prefetch_forward(&spinorfield->BGQ_CONSPTR[d][beginj]);
	}

	for (size_t j = beginj; j < endj; j+=1) {
#ifndef NDEBUG
		size_t offset = bgq_pointer2offset(spinorfield, &g_bgq_sec_recv[d][j]);
		ucoord index = bgq_offset2index(offset);
		ucoord ic = bgq_index2collapsed(isOdd, index, -1);
		ucoord t1 = bgq_collapsed2t(isOdd, ic, 0);
		ucoord t2 = bgq_collapsed2t(isOdd, ic, 1);
		ucoord x = bgq_collapsed2x(isOdd, ic);
		ucoord y = bgq_collapsed2y(isOdd, ic);
		ucoord z = bgq_collapsed2z(isOdd, ic);
		assert(bgq_offset2ddst(offset)==d);
#endif

		bgq_weyl_vec *weyladdr_src = &g_bgq_sec_recv[d][j];
		bgq_ptrnext_prefetch(&spinorfield->BGQ_CONSPTR[d][j]);
		bgq_weyl_vec *weyladdr_dst = spinorfield->BGQ_CONSPTR[isOdd][d][j];

		bgq_su3_weyl_decl(weyl);
		bgq_su3_weylnext_prefetch(weyladdr_src);
		bgq_su3_weyl_load(weyl, weyladdr_src);
				bgq_weylqpx_expect(weyl, t1, t2, x, y, z, d, false);

		bgq_su3_weyl_store(weyladdr_dst, weyl);
	}
}


static inline void bgq_HoppingMatrix_worker_datamove_recvtup_unvectorized(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t begini, size_t endi, bool noprefetchstream) {
	assert(COMM_T);
	assert(BGQ_UNVECTORIZE);

	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_sec_recv[TUP][begini]);
		bgq_prefetch_forward(&g_bgq_sec_temp_tdown[2*begini]);
	}

	for (size_t i = begini; i < endi; i+=1) {
		bgq_weyl_vec *weyladdr_right = &g_bgq_sec_recv[TUP][i];

		bgq_su3_weyl_decl(weyl_right);
		bgq_su3_weylnext_prefetch(weyladdr_right);
		bgq_su3_weyl_load(weyl_right, weyladdr_right);

		{
			size_t j = 2*i;
#ifndef NDEBUG
			size_t offset_left = bgq_pointer2offset(spinorfield, &g_bgq_sec_temp_tdown[j]);
			ucoord index_left = bgq_offset2index(offset_left);
			size_t offset_right = bgq_pointer2offset(spinorfield, weyladdr_right);
			ucoord index_right = bgq_offset2index(offset_right);
			ucoord ic_left = bgq_index2collapsed(isOdd, index_left, -1);
			ucoord ic_right = bgq_index2collapsed(isOdd, index_right, 0);
			assert(ic_left == ic_right);
			ucoord ic = ic_left;
			ucoord t1 = bgq_collapsed2t(isOdd, ic_left, 0);
			ucoord t2 = bgq_collapsed2t(isOdd, ic_right, 1);
			ucoord x = bgq_collapsed2x(isOdd, ic);
			ucoord y = bgq_collapsed2y(isOdd, ic);
			ucoord z = bgq_collapsed2z(isOdd, ic);
			bgq_direction d1 = bgq_offset2ddst(offset_left);
			bgq_direction d2 = bgq_offset2ddst(offset_right);
			assert(d1 == d2);
			bgq_direction d = d1;
#endif
					bgq_weylqpxk_expect(weyl_right, 0, t2, x, y, z, d, false);

			//TODO: Is reading just the 16 used bytes faster?
			bgq_weyl_vec *weyladdr_left = &g_bgq_sec_temp_tdown[j];
			bgq_ptrnext_prefetch(&spinorfield->BGQ_CONSPTR[TUP][j]);
			bgq_weyl_vec *weyladdr_dst = spinorfield->BGQ_CONSPTR[isOdd][TUP][j];

			bgq_su3_weyl_decl(weyl_left);
			bgq_su3_weylnext_prefetch(weyladdr_left);
			bgq_su3_weyl_load(weyl_left, weyladdr_left);
					bgq_weylqpxk_expect(weyl_left, 1, t1, x, y, z, d1, false);

			bgq_su3_weyl_decl(weyl);
			bgq_su3_weyl_merge2(weyl, weyl_left, weyl_right);
					bgq_weylqpx_expect(weyl, t1, t2, x, y, z, d, false);

#ifndef NDEBUG
			size_t offset = bgq_pointer2offset(spinorfield, weyladdr_dst);
			ucoord index = bgq_offset2index(offset);
			bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
			assert(sec==sec_body || sec==sec_surface);
			ucoord ic_check = g_bgq_index2collapsed[isOdd][index];
			assert(ic == ic_check);
#endif
			bgq_su3_weyl_store(weyladdr_dst, weyl);
		}

		{
			size_t j = 2*i + 1;
#ifndef NDEBUG
			size_t offset_left = bgq_pointer2offset(spinorfield, &g_bgq_sec_temp_tdown[j]);
			ucoord index_left = bgq_offset2index(offset_left);
			size_t offset_right = bgq_pointer2offset(spinorfield, weyladdr_right);
			ucoord index_right = bgq_offset2index(offset_right);
			ucoord ic_left =  bgq_index2collapsed(isOdd, index_left, -1);
			ucoord ic_right = bgq_index2collapsed(isOdd, index_right, 1);
			assert(ic_left == ic_right);
			ucoord ic = ic_left;
			ucoord t1 = bgq_collapsed2t(isOdd, ic_left, 0);
			ucoord t2 = bgq_collapsed2t(isOdd, ic_right, 1);
			ucoord x = bgq_collapsed2x(isOdd, ic);
			ucoord y = bgq_collapsed2y(isOdd, ic);
			ucoord z = bgq_collapsed2z(isOdd, ic);
			bgq_direction d1 = bgq_offset2ddst(offset_left);
			bgq_direction d2 = bgq_offset2ddst(offset_right);
			assert(d1==d2);
			bgq_direction d = d1;
#endif
					bgq_weylqpxk_expect(weyl_right, 1, t2, x, y, z, d, false);

			//TODO: Is reading just the 16 used bytes faster?
			bgq_weyl_vec *weyladdr_left = &g_bgq_sec_temp_tdown[j];
			bgq_ptrnext_prefetch(&spinorfield->BGQ_CONSPTR[TUP][j]);
			bgq_weyl_vec *weyladdr_dst = spinorfield->BGQ_CONSPTR[isOdd][TUP][j];

			bgq_su3_weyl_decl(weyl_left);
			bgq_su3_weyl_load(weyl_left, weyladdr_left);
			bgq_su3_weylnext_prefetch(weyladdr_left);
					bgq_weylqpxk_expect(weyl_left, 1, t1, x, y, z, d1, false);

			bgq_su3_weyl_decl(weyl);
			bgq_su3_weyl_rmerge(weyl, weyl_left, weyl_right);
					bgq_weylqpx_expect(weyl, t1, t2, x, y, z, d, false);

#ifndef NDEBUG
			size_t offset = bgq_pointer2offset(spinorfield, weyladdr_dst);
			ucoord index = bgq_offset2index(offset);
			bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
			assert(sec==sec_body || sec==sec_surface);
			ucoord ic_check = g_bgq_index2collapsed[isOdd][index];
			assert(ic == ic_check);
#endif
			bgq_su3_weyl_store(weyladdr_dst, weyl);
		}
	}
}


static inline void bgq_HoppingMatrix_worker_datamove_recvtdown_unvectorized(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t begini, size_t endi, bool noprefetchstream) {
	assert(COMM_T);
	assert(BGQ_UNVECTORIZE);

	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_sec_recv[TDOWN][begini]);
		bgq_prefetch_forward(&g_bgq_sec_temp_tup[2*begini]);
	}


	for (size_t i = begini; i < endi; i+=1) {
		bgq_weyl_vec *weyladdr_left = &g_bgq_sec_recv[TDOWN][i];

		bgq_su3_weyl_decl(weyl_left);
		bgq_su3_weylnext_prefetch(weyladdr_left);
		bgq_su3_weyl_load(weyl_left, weyladdr_left);

		{
			size_t j = 2*i;
#ifndef NDEBUG
			size_t offset_left = bgq_pointer2offset(spinorfield, weyladdr_left);
			ucoord index_left = bgq_offset2index(offset_left);
			size_t offset_right = bgq_pointer2offset(spinorfield, &g_bgq_sec_temp_tup[j]);
			ucoord index_right = bgq_offset2index(offset_right);
			ucoord ic_left = bgq_index2collapsed(isOdd, index_left, 0);
			ucoord ic_right = bgq_index2collapsed(isOdd, index_right, -1);
			assert(ic_left == ic_right);
			ucoord ic = ic_left;
			ucoord t1 = bgq_collapsed2t(isOdd, ic_left, 0);
			ucoord t2 = bgq_collapsed2t(isOdd, ic_right, 1);
			ucoord x = bgq_collapsed2x(isOdd, ic);
			ucoord y = bgq_collapsed2y(isOdd, ic);
			ucoord z = bgq_collapsed2z(isOdd, ic);
			bgq_direction d1 = bgq_offset2ddst(offset_left);
			bgq_direction d2 = bgq_offset2ddst(offset_right);
			assert(d1 == d2);
			bgq_direction d = d1;
#endif
					bgq_weylqpxk_expect(weyl_left, 0, t1, x, y, z, d, false);

			//TODO: Is reading just the 16 used bytes faster?
			bgq_weyl_vec *weyladdr_right = &g_bgq_sec_temp_tup[j];
			bgq_ptrnext_prefetch(&spinorfield->BGQ_CONSPTR[TDOWN][j]);
			bgq_weyl_vec *weyladdr_dst = spinorfield->BGQ_CONSPTR[isOdd][TDOWN][j];

			bgq_su3_weyl_decl(weyl_right);
			bgq_su3_weylnext_prefetch(weyladdr_right);
			bgq_su3_weyl_load(weyl_right, weyladdr_right);
					bgq_weylqpxk_expect(weyl_right, 0, t2, x, y, z, d1, false);

			bgq_su3_weyl_decl(weyl);
			bgq_su3_weyl_lmerge(weyl, weyl_left, weyl_right);
					bgq_weylqpx_expect(weyl, t1, t2, x, y, z, d, false);

			bgq_su3_weyl_store(weyladdr_dst, weyl);
		}

		{
			size_t j = 2*i + 1;
#ifndef NDEBUG
			size_t offset_left = bgq_pointer2offset(spinorfield, weyladdr_left);
			ucoord index_left = bgq_offset2index(offset_left);
			size_t offset_right = bgq_pointer2offset(spinorfield, &g_bgq_sec_temp_tup[j]);
			ucoord index_right = bgq_offset2index(offset_right);
			ucoord ic_left = bgq_index2collapsed(isOdd, index_left, 1);
			ucoord ic_right = bgq_index2collapsed(isOdd, index_right, -1);
			assert(ic_left == ic_right);
			ucoord ic = ic_left;
			ucoord t1 = bgq_collapsed2t(isOdd, ic_left, 0);
			ucoord t2 = bgq_collapsed2t(isOdd, ic_right, 1);
			ucoord x = bgq_collapsed2x(isOdd, ic);
			ucoord y = bgq_collapsed2y(isOdd, ic);
			ucoord z = bgq_collapsed2z(isOdd, ic);
			bgq_direction d1 = bgq_offset2ddst(offset_left);
			bgq_direction d2 = bgq_offset2ddst(offset_right);
			assert(d1 == d2);
			bgq_direction d = d1;
#endif
					bgq_weylqpxk_expect(weyl_left, 1, t1, x, y, z, d, false);

			//TODO: Is reading just the 16 used bytes faster?
			bgq_weyl_vec *weyladdr_right = &g_bgq_sec_temp_tup[j];
			bgq_ptrnext_prefetch(&spinorfield->BGQ_CONSPTR[TDOWN][j]);
			bgq_weyl_vec *weyladdr_dst = spinorfield->BGQ_CONSPTR[isOdd][TDOWN][j];

			bgq_su3_weyl_decl(weyl_right);
			bgq_su3_weylnext_prefetch(weyladdr_right);
			bgq_su3_weyl_load(weyl_right, weyladdr_right);
					bgq_weylqpxk_expect(weyl_right, 0, t2, x, y, z, d1, false);

			bgq_su3_weyl_decl(weyl);
			bgq_su3_weyl_merge2(weyl, weyl_left, weyl_right);
					bgq_weylqpx_expect(weyl, t1, t2, x, y, z, d, false);

			bgq_su3_weyl_store(weyladdr_dst, weyl);
		}
	}
}


static inline void bgq_HoppingMatrix_worker_datamove_recvtup(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t beginj, size_t endj, bool noprefetchstream) {
	assert(COMM_T);
	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_sec_send[TDOWN][beginj]);
		bgq_prefetch_forward(&g_bgq_sec_recv[TUP][beginj]);
		bgq_prefetch_forward(&spinorfield->BGQ_CONSPTR[TUP][beginj]);
	}

	for (size_t j = beginj; j < endj; j+=1) {
#ifndef NDEBUG
		size_t offset_left = bgq_pointer2offset(spinorfield, &g_bgq_sec_send[TDOWN][j]);
		ucoord index_left = bgq_offset2index(offset_left);
		size_t offset_right =bgq_pointer2offset(spinorfield, &g_bgq_sec_recv[TUP][j]);
		ucoord index_right = bgq_offset2index(offset_right);
		ucoord ic_left = bgq_index2collapsed(isOdd, index_left, -1);
		ucoord ic_right = bgq_index2collapsed(isOdd, index_right, -1);
		assert(ic_left == ic_right);
		ucoord ic = ic_left;
		ucoord t1 = bgq_collapsed2t(isOdd, ic_left, 0);
		ucoord t2 = bgq_collapsed2t(isOdd, ic_right, 1);
		ucoord x = bgq_collapsed2x(isOdd, ic);
		ucoord y = bgq_collapsed2y(isOdd, ic);
		ucoord z = bgq_collapsed2z(isOdd, ic);
		bgq_direction d1 = bgq_offset2ddst(offset_left);
		bgq_direction d2 = bgq_offset2ddst(offset_right);
		assert(d1==d2);
		bgq_direction d = d1;
#endif

		//TODO: Is reading just the 16 used bytes faster?
		bgq_weyl_vec *weyladdr_left = &g_bgq_sec_send[TDOWN][j];
		bgq_weyl_vec *weyladdr_right = &g_bgq_sec_recv[TUP][j];
		bgq_ptrnext_prefetch(&spinorfield->BGQ_CONSPTR[TUP][j]);
		bgq_weyl_vec *weyladdr_dst = spinorfield->BGQ_CONSPTR[isOdd][TUP][j];

		bgq_su3_weyl_decl(weyl_left);
		bgq_su3_weylnext_prefetch(weyladdr_left);
		bgq_su3_weyl_load(weyl_left, weyladdr_left);
		assert(bgq_cmplxval2(weyl_left_v0_c0)!=0); // for valgrind
				bgq_weylqpxk_expect(weyl_left, 1, t1, x, y, z, d1, false);

		bgq_su3_weyl_decl(weyl_right);
		bgq_su3_weylnext_prefetch(weyladdr_right);
		bgq_su3_weyl_load(weyl_right, weyladdr_right);
		assert(bgq_cmplxval1(weyl_left_v0_c0)!=0); // for valgrind
				bgq_weylqpxk_expect(weyl_right, 0, t2, x, y, z, d2, false);

		bgq_su3_weyl_decl(weyl);
		bgq_su3_weyl_merge2(weyl, weyl_left, weyl_right);
				bgq_weylqpx_expect(weyl, t1, t2, x, y, z, d, false);

#ifndef NDEBUG
		size_t offset = bgq_pointer2offset(spinorfield, weyladdr_dst);
		ucoord index = bgq_offset2index(offset);
		bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
		assert(sec==sec_body || sec==sec_surface);
		ucoord ic_check = g_bgq_index2collapsed[isOdd][index];
		assert(ic == ic_check);
#endif
		bgq_su3_weyl_store(weyladdr_dst, weyl);
	}
}


static inline void bgq_HoppingMatrix_worker_datamove_recvtdown(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t beginj, size_t endj, bool noprefetchstream) {
	assert(COMM_T);
	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_sec_recv[TDOWN][beginj]);
		bgq_prefetch_forward(&g_bgq_sec_send[TUP][beginj]);
		bgq_prefetch_forward(&spinorfield->BGQ_CONSPTR[TDOWN][beginj]);
	}

	for (size_t j = beginj; j < endj; j+=1) {
#ifndef NDEBUG
		size_t offset_left = bgq_pointer2offset(spinorfield, &g_bgq_sec_recv[TDOWN][j]);
		ucoord index_left = bgq_offset2index(offset_left);
		size_t offset_right = bgq_pointer2offset(spinorfield, &g_bgq_sec_send[TUP][j]);
		ucoord index_right = bgq_offset2index(offset_right);
		ucoord ic_left = bgq_index2collapsed(isOdd, index_left, -1);
		ucoord ic_right = bgq_index2collapsed(isOdd, index_right, -1);
		assert(ic_left == ic_right);
		ucoord ic = ic_left;
		ucoord t1 = bgq_collapsed2t(isOdd, ic_left, 0);
		ucoord t2 = bgq_collapsed2t(isOdd, ic_right, 1);
		ucoord x = bgq_collapsed2x(isOdd, ic);
		ucoord y = bgq_collapsed2y(isOdd, ic);
		ucoord z = bgq_collapsed2z(isOdd, ic);
		bgq_direction d1 = bgq_offset2ddst(offset_left);
		bgq_direction d2 = bgq_offset2ddst(offset_right);
		assert(d1==d2);
		bgq_direction d = d1;
#endif

		bgq_weyl_vec *weyladdr_left = &g_bgq_sec_recv[TDOWN][j];
		bgq_weyl_vec *weyladdr_right = &g_bgq_sec_send[TUP][j];
		bgq_ptrnext_prefetch(&spinorfield->BGQ_CONSPTR[TDOWN][j]);
		bgq_weyl_vec *weyladdr_dst = spinorfield->BGQ_CONSPTR[isOdd][TDOWN][j];

		bgq_su3_weyl_decl(weyl_left);
		bgq_su3_weylnext_prefetch(weyladdr_left);
		bgq_su3_weyl_load(weyl_left, weyladdr_left);
				bgq_weylqpxk_expect(weyl_left, 1, t1, x, y, z, d1, false);

		bgq_su3_weyl_decl(weyl_right);
		bgq_su3_weylnext_prefetch(weyladdr_right);
		bgq_su3_weyl_load(weyl_right, weyladdr_right);
				bgq_weylqpxk_expect(weyl_right, 0, t2, x, y, z, d2, false);

		bgq_su3_weyl_decl(weyl);
		bgq_su3_weyl_merge2(weyl, weyl_left, weyl_right);
				bgq_weylqpx_expect(weyl, t1, t2, x, y, z, d, false);

#ifndef NDEBUG
		size_t offset = bgq_pointer2offset(spinorfield, weyladdr_dst);
		ucoord index = bgq_offset2index(offset);
		bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
		assert(sec==sec_body || sec==sec_surface);
		ucoord ic_check = g_bgq_index2collapsed[isOdd][index];
		assert(ic == ic_check);
#endif
		bgq_su3_weyl_store(weyladdr_dst, weyl);
	}
}


void bgq_HoppingMatrix_worker_datamove(void *arg_untyped, size_t tid, size_t threads) {
	bgq_work_datamove *arg = arg_untyped;
	bgq_weylfield_controlblock *spinorfield = arg->spinorfield;
	bgq_hmflags opts = arg->opts;
	bool isOdd = spinorfield->isOdd;
	bool noprefetchstream = opts & hm_noprefetchstream;

	assert(COMM_T);
#if BGQ_UNVECTORIZE
	const size_t workload_recv_tup = COMM_T * LOCAL_HALO_T/(PHYSICAL_LP*PHYSICAL_LK);
#else
	const size_t workload_recv_tup =  COMM_T * LOCAL_HALO_T/PHYSICAL_LP;
#endif
	const size_t workload_recv_tdown = workload_recv_tup;
	const size_t workload_recv = 2*PHYSICAL_HALO_X + 2*PHYSICAL_HALO_Y + 2*PHYSICAL_HALO_Z;
	const size_t workload = workload_recv + 2*workload_recv_tup + 2*workload_recv_tdown;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (size_t i = begin; i<end; ) {
		WORKLOAD_DECL(i, workload);

		if (WORKLOAD_SPLIT(workload_recv)) {
			size_t seclength;
			bgq_dimension dim;
			bool isDown = WORKLOAD_CHUNK(2);
			if (WORKLOAD_SPLIT(PHYSICAL_HALO_X)) {
				seclength = PHYSICAL_HALO_X;
				dim = DIM_X;
			} else if (WORKLOAD_SPLIT(PHYSICAL_HALO_Y)) {
				seclength = PHYSICAL_HALO_Y;
				dim = DIM_Y;
			} else if (WORKLOAD_SPLIT(PHYSICAL_HALO_Z)) {
				seclength = PHYSICAL_HALO_Z;
				dim = DIM_Z;
			} else {
				UNREACHABLE
			}
			size_t beginj = WORKLOAD_PARAM(seclength);
			size_t endj = min_sizet(seclength,beginj+threadload);
			bgq_direction d = bgq_direction_compose(dim, isDown);
			bgq_HoppingMatrix_worker_datamove_recvxyz(spinorfield, d, isOdd, beginj, endj, noprefetchstream);
			i += (endj - beginj);
		} else if (WORKLOAD_SPLIT(2*workload_recv_tup)) {
			// Count an T-iteration twice; better have few underloaded threads (so remaining SMT-threads have some more resources) then few overloaded threads (so the master thread has to wait for them)
			size_t twobeginj = WORKLOAD_PARAM(2*workload_recv_tup);
			size_t twoendj = min_sizet(2*workload_recv_tup, twobeginj+threadload);
			size_t beginj = twobeginj / 2;
			size_t endj = twoendj / 2;
			if (BGQ_UNVECTORIZE)
				bgq_HoppingMatrix_worker_datamove_recvtup_unvectorized(spinorfield, isOdd, beginj, endj, noprefetchstream);
			else
				bgq_HoppingMatrix_worker_datamove_recvtup(spinorfield, isOdd, beginj, endj, noprefetchstream);
			i += (twoendj - twobeginj);
		} else if (WORKLOAD_SPLIT(2*workload_recv_tdown)) {
			size_t twobeginj = WORKLOAD_PARAM(2*workload_recv_tdown);
			size_t twoendj = min_sizet(2*workload_recv_tdown, twobeginj+threadload);
			size_t beginj = twobeginj / 2;
			size_t endj = twoendj / 2;
			if (BGQ_UNVECTORIZE)
				bgq_HoppingMatrix_worker_datamove_recvtdown_unvectorized(spinorfield, isOdd, beginj, endj, noprefetchstream);
			else
				bgq_HoppingMatrix_worker_datamove_recvtdown(spinorfield, isOdd, beginj, endj, noprefetchstream);
			i += (twoendj - twobeginj);
		} else {
			UNREACHABLE
		}
		WORKLOAD_CHECK
	}
}


static inline void bgq_HoppingMatrix_worker_datamovet_recvtup(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t beginj, size_t endj, bool noprefetchstream) {
	assert(!COMM_T);
	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_sec_temp_tdown[beginj]);
		bgq_prefetch_forward(&spinorfield->BGQ_CONSPTR[TUP][beginj]);
	}

	for (size_t j = beginj; j < endj; j+=1) {
#ifndef NDEBUG
		size_t offset = bgq_pointer2offset(spinorfield, &g_bgq_sec_temp_tdown[j]);
		ucoord index = bgq_offset2index(offset);
		ucoord ic = bgq_index2collapsed(isOdd, index, -1);
		ucoord t1 = bgq_collapsed2t(isOdd, ic, 0);
		ucoord t2 = bgq_collapsed2t(isOdd, ic, 1);
		ucoord x = bgq_collapsed2x(isOdd, ic);
		ucoord y = bgq_collapsed2y(isOdd, ic);
		ucoord z = bgq_collapsed2z(isOdd, ic);
		bgq_direction d = bgq_offset2ddst(offset);
#endif

		bgq_weyl_vec *weyladdr = &g_bgq_sec_temp_tdown[j];
		bgq_ptrnext_prefetch(&spinorfield->BGQ_CONSPTR[TUP][j]);
		bgq_weyl_vec *weyladdr_dst = spinorfield->BGQ_CONSPTR[isOdd][TUP][j];
		assert(weyladdr_dst);

		bgq_su3_weyl_decl(weyl_before);
		bgq_su3_weylnext_prefetch(weyladdr);
		bgq_su3_weyl_load(weyl_before, weyladdr);
				bgq_weylqpx_expect(weyl_before, t2, t1, x, y, z, d, false);

		bgq_su3_weyl_decl(weyl_after);
		bgq_su3_weyl_merge2(weyl_after, weyl_before, weyl_before); // Just switch first and second complex
				bgq_weylqpx_expect(weyl_after, t1, t2, x, y, z, d, false);

		bgq_su3_weyl_store(weyladdr_dst, weyl_after);
	}
}


static inline void bgq_HoppingMatrix_worker_datamovet_recvtdown(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t beginj, size_t endj, bool noprefetchstream) {
	assert(!COMM_T);
	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_sec_temp_tup[beginj]);
		bgq_prefetch_forward(&spinorfield->BGQ_CONSPTR[TDOWN][beginj]);
	}

	for (size_t j = beginj; j < endj; j+=1) {
#ifndef NDEBUG
		size_t offset = bgq_pointer2offset(spinorfield, &g_bgq_sec_temp_tup[j]);
		ucoord index = bgq_offset2index(offset);
		ucoord ic = g_bgq_index2collapsed[isOdd][index];
		ucoord t1 = bgq_collapsed2t(isOdd, ic, 0);
		ucoord t2 = bgq_collapsed2t(isOdd, ic, 1);
		ucoord x = bgq_collapsed2x(isOdd, ic);
		ucoord y = bgq_collapsed2y(isOdd, ic);
		ucoord z = bgq_collapsed2z(isOdd, ic);
		bgq_direction d = bgq_offset2ddst(offset);
#endif

		bgq_weyl_vec *weyladdr = &g_bgq_sec_temp_tup[j];
		bgq_ptrnext_prefetch(&spinorfield->BGQ_CONSPTR[TDOWN][j]);
		bgq_weyl_vec *weyladdr_dst = spinorfield->BGQ_CONSPTR[isOdd][TDOWN][j];
		assert(weyladdr_dst);

		bgq_su3_weyl_decl(weyl_before);
		bgq_su3_weylnext_prefetch(weyladdr);
		bgq_su3_weyl_load(weyl_before, weyladdr);
				bgq_weylqpx_expect(weyl_before, t2, t1, x, y, z, d, false);

		bgq_su3_weyl_decl(weyl_after);
		bgq_su3_weyl_merge2(weyl_after, weyl_before, weyl_before);
				bgq_weylqpx_expect(weyl_after, t1, t2, x, y, z, d, false);

		bgq_su3_weyl_store(weyladdr_dst, weyl_after);
	}
}


void bgq_HoppingMatrix_datamovet_worker(void *arg_untyped, size_t tid, size_t threads) {
	assert(!COMM_T);

	bgq_work_datamove *arg = arg_untyped;
	bgq_weylfield_controlblock *spinorfield = arg->spinorfield;
	bgq_hmflags opts = arg->opts;
	bool isOdd = spinorfield->isOdd;
	bool noprefetchstream = opts & hm_noprefetchstream;

	const size_t workload_recv_tup = LOCAL_HALO_T/PHYSICAL_LP;
	const size_t workload_recv_tdown = workload_recv_tup;
	const size_t workload = workload_recv_tup + workload_recv_tdown;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (size_t i = begin; i<end; ) {
		WORKLOAD_DECL(i,workload);

		if (WORKLOAD_SPLIT(workload_recv_tup)) {
			size_t beginj = WORKLOAD_PARAM(workload_recv_tup);
			size_t endj = min_sizet(workload_recv_tup, beginj+threadload);
			bgq_HoppingMatrix_worker_datamovet_recvtup(spinorfield, isOdd, beginj, endj, noprefetchstream);
			i += (endj - beginj);
		} else if (WORKLOAD_SPLIT(workload_recv_tdown)) {
			size_t beginj = WORKLOAD_PARAM(workload_recv_tdown);
			size_t endj = min_sizet(workload_recv_tdown,beginj+threadload);
			bgq_HoppingMatrix_worker_datamovet_recvtdown(spinorfield, isOdd, beginj, endj, noprefetchstream);
			i += (endj - beginj);
		} else {
			UNREACHABLE
		}

		WORKLOAD_CHECK
	}
}


void bgq_copyFromLegacy_worker(void *arg_untyped, size_t tid, size_t threads) {
	bgq_copyFromLegacy_workload *arg = arg_untyped;
	bool isOdd = arg->isOdd;
	spinor *sourcefield = arg->source;
	bgq_weylfield_controlblock *targetfield = arg->target;
	assert(targetfield->BGQ_SEC_FULLLAYOUT);

	const size_t workload = VOLUME/2;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);

	bgq_prefetch_forward(&sourcefield[begin]);

	for (int eosub = begin; eosub<end; eosub+=1) {
		spinor *srcarddr = &sourcefield[eosub];

		ucoord ic = bgq_eosub2collapsed(isOdd, eosub);
		ucoord k = bgq_eosub2k(isOdd, eosub);
#ifndef NDEBUG
		ucoord t = bgq_collapsed2t(isOdd, ic, k);
		ucoord x = bgq_collapsed2x(isOdd, ic);
		ucoord y = bgq_collapsed2y(isOdd, ic);
		ucoord z = bgq_collapsed2z(isOdd, ic);
#endif
		bgq_su3_spinor_decl(spinor);
		bgq_su3_spinornext_prefetch_double(srcarddr);
		bgq_ld2a_double(spinor_v0_c0, 0, srcarddr);
		bgq_qvlfcduxa(spinor_v0_c1, srcarddr, 16);
		bgq_qvlfcduxa(spinor_v0_c2, srcarddr, 16);
		bgq_qvlfcduxa(spinor_v1_c0, srcarddr, 16);
		bgq_qvlfcduxa(spinor_v1_c1, srcarddr, 16);
		bgq_qvlfcduxa(spinor_v1_c2, srcarddr, 16);
		bgq_qvlfcduxa(spinor_v2_c0, srcarddr, 16);
		bgq_qvlfcduxa(spinor_v2_c1, srcarddr, 16);
		bgq_qvlfcduxa(spinor_v2_c2, srcarddr, 16);
		bgq_qvlfcduxa(spinor_v3_c0, srcarddr, 16);
		bgq_qvlfcduxa(spinor_v3_c1, srcarddr, 16);
		bgq_qvlfcduxa(spinor_v3_c2, srcarddr, 16);
				bgq_spinorqpxk_expect(spinor,0,t,x,y,z);

		bgq_spinor_vec *targetaddr = &targetfield->BGQ_SEC_FULLLAYOUT[ic];
		targetaddr = (bgq_spinor_vec*)((char*)targetaddr + k*PRECISION_COMPLEX_SIZEOF);
		bgq_st2a(spinor_v0_c0,0,targetaddr);
		bgq_qvstfcuxa(spinor_v0_c1,targetaddr,2*PRECISION_COMPLEX_SIZEOF); // i.e. skip 16 bytes, write 16 bytes
		bgq_qvstfcuxa(spinor_v0_c2,targetaddr,2*PRECISION_COMPLEX_SIZEOF);
		bgq_qvstfcuxa(spinor_v1_c0,targetaddr,2*PRECISION_COMPLEX_SIZEOF);
		bgq_qvstfcuxa(spinor_v1_c1,targetaddr,2*PRECISION_COMPLEX_SIZEOF);
		bgq_qvstfcuxa(spinor_v1_c2,targetaddr,2*PRECISION_COMPLEX_SIZEOF);
		bgq_qvstfcuxa(spinor_v2_c0,targetaddr,2*PRECISION_COMPLEX_SIZEOF);
		bgq_qvstfcuxa(spinor_v2_c1,targetaddr,2*PRECISION_COMPLEX_SIZEOF);
		bgq_qvstfcuxa(spinor_v2_c2,targetaddr,2*PRECISION_COMPLEX_SIZEOF);
		bgq_qvstfcuxa(spinor_v3_c0,targetaddr,2*PRECISION_COMPLEX_SIZEOF);
		bgq_qvstfcuxa(spinor_v3_c1,targetaddr,2*PRECISION_COMPLEX_SIZEOF);
		bgq_qvstfcuxa(spinor_v3_c2,targetaddr,2*PRECISION_COMPLEX_SIZEOF);
				bgq_spinorveck_expect(targetfield->BGQ_SEC_FULLLAYOUT[ic], k, t, x, y, z);
	}
}


static inline void bgq_spinorfield_rewrite_worker(void *arg_untyped, size_t tid, size_t threads, bool weyllayout, bool sloppy, bool mul, bool isLegacy) {
	bgq_spinorfield_rewrite_work *arg = arg_untyped;
	bgq_weylfield_controlblock *field  = arg->field;
	bool isOdd = arg->isOdd;

	assert(field->BGQ_SEC_FULLLAYOUT);
	assert(field->BGQ_HAS_FULLLAYOUT);
	const size_t workload = PHYSICAL_VOLUME;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);

	bgq_spinorfield_streamSpinor(field, isOdd, begin, weyllayout, sloppy, mul, isLegacy);
	for (ucoord ic = begin; ic < end; ic+=1) {
#ifndef NDEBUG
		ucoord ih = bgq_collapsed2halfvolume(field->isOdd, ic);
		ucoord t1 = bgq_halfvolume2t1(field->isOdd, ih);
		ucoord t2 = bgq_halfvolume2t2(field->isOdd, ih);
		ucoord tv = bgq_halfvolume2tv(ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);
#endif

		bgq_su3_spinor_decl(spinor);
		bgq_spinorfield_readSpinor(&spinor, field, isOdd, ic, weyllayout, sloppy, mul, isLegacy);
		bgq_spinorfield_prefetchNextSpinor(field, isOdd, ic, weyllayout, sloppy, mul, isLegacy);

		bgq_spinor_vec *fulladdr = &field->BGQ_SEC_FULLLAYOUT[ic];
		bgq_su3_spinor_store(fulladdr, spinor);
	}
}


BGQ_SPINORFIELD_GENWORKER(bgq_spinorfield_rewrite_worker)
