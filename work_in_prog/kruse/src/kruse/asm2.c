
#include "bgq_HoppingMatrix.h"

#include "bgq_field.h"
#include "bgq_spinorfield.h"
#include "bgq_qpx.h"
#include "bgq_dispatch.h"
#include "bgq_comm.h"
#include "bgq_workers.h"
#include "bgq_gaugefield.h"

#include "../boundary.h"
#include "../update_backward_gauge.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>


static inline void bgq_site_rmul_plain_sub(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), bgq_params(qr), ucoord ic) {
	bgq_su3_spinor_decl(result);
	bgq_su3_spinor_decl(rmul1);

	// Unfortunately, there is no fused-multiply-sub for this
	bgq_mul(rmul1_v0_c0, qr, spinor1_v0_c0);
	bgq_mul(rmul1_v0_c1, qr, spinor1_v0_c1);
	bgq_mul(rmul1_v0_c2, qr, spinor1_v0_c2);
	bgq_mul(rmul1_v1_c0, qr, spinor1_v1_c0);
	bgq_mul(rmul1_v1_c1, qr, spinor1_v1_c1);
	bgq_mul(rmul1_v1_c2, qr, spinor1_v1_c2);
	bgq_mul(rmul1_v2_c0, qr, spinor1_v2_c0);
	bgq_mul(rmul1_v2_c1, qr, spinor1_v2_c1);
	bgq_mul(rmul1_v2_c2, qr, spinor1_v2_c2);
	bgq_mul(rmul1_v3_c0, qr, spinor1_v3_c0);
	bgq_mul(rmul1_v3_c1, qr, spinor1_v3_c1);
	bgq_mul(rmul1_v3_c2, qr, spinor1_v3_c2);

	bgq_su3_vsub(result_v0, rmul1_v0, spinor2_v0);
	bgq_su3_vsub(result_v1, rmul1_v1, spinor2_v1);
	bgq_su3_vsub(result_v2, rmul1_v2, spinor2_v2);
	bgq_su3_vsub(result_v3, rmul1_v3, spinor2_v3);

	bgq_su3_spinor_mov(*target, result);
}


/*
 * bgq_operator.inc.c
 *
 *  Created on: Nov 9, 2012
 *      Author: meinersbur
 */

#include "bgq_spinorfield.h"
#include "bgq_dispatch.h"
#include "bgq_HoppingMatrix.h"
#include "bgq_qpx.h"
#include "bgq_utils.h"

#include <stdbool.h>


#ifndef OPERATOR_NAME
// Make this .inc.c file compilable as stand-alone
#define OPERATOR_NAME bgq_spinorfield_rmul_plain_sub_raw
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_rmul_plain_sub
#define OPERATOR_EXTRAPARMS bgq_params(qr)
#define OPERATOR_EXTRAARGS bgq_vars(qr)


#endif


#ifndef OPERATOR_EXTRAPARMS
#define OPERATOR_EXTRAPARMS
#endif
#ifndef OPERATOR_EXTRAARGS
#define OPERATOR_EXTRAARGS
#endif

#define OPERATOR_EXTRAPARMLIST IF_EMPTY_INPAREN((),(,),OPERATOR_EXTRAPARMS) OPERATOR_EXTRAPARMS
#define OPERATOR_EXTRAARGLIST IF_EMPTY_INPAREN((),(,),OPERATOR_EXTRAARGS) OPERATOR_EXTRAARGS
#define OPERATOR_EXTRA_ASSIGN(varname) (varname) = ((arg->extra).varname);

#if OPERATOR_ARGFIELDS==0
#define IFNOARG(...) __VA_ARGS__
#else
#define IFNOARG(...)
#endif

#if OPERATOR_ARGFIELDS>=1
#define IF1ARG(...) __VA_ARGS__
#else
#define IF1ARG(...)
#endif

#if OPERATOR_ARGFIELDS>=2
#define IF2ARG(...) __VA_ARGS__
#else
#define IF2ARG(...)
#endif

#if ISEMPTY(OPERATOR_EXTRAPARMS)
#define IFEXTRA(...)
#else
#define IFEXTRA(...) __VA_ARGS__
#endif

#define OPERATOR_LAYOUTNAME_ly_full_double Fulllayout
#define OPERATOR_LAYOUTNAME_ly_weyl_double Weyllayout
#define OPERATOR_LAYOUTNAME_ly_full_float FulllayoutSloppy
#define OPERATOR_LAYOUTNAME_ly_weyl_float WeyllayoutSloppy
#define OPERATOR_LAYOUTNAME_0 OPERATOR_LAYOUTNAME_ly_full_double
#define OPERATOR_LAYOUTNAME_1 OPERATOR_LAYOUTNAME_ly_weyl_double
#define OPERATOR_LAYOUTNAME_2 OPERATOR_LAYOUTNAME_ly_full_float
#define OPERATOR_LAYOUTNAME_3 OPERATOR_LAYOUTNAME_ly_weyl_float

#define OPERATOR_PRECISION_ISSLOPPY_double false
#define OPERATOR_PRECISION_ISSLOPPY_float true


#if !ISEMPTY(OPERATOR_EXTRAPARMS)
typedef struct {
	PREPROCESSOR_FOREACH( ,;,;, ,IDENTITY,OPERATOR_EXTRAPARMS)
} NAME2(OPERATOR_NAME,extraparm_t);
#endif

typedef struct {
	tristate isOdd;
	bgq_weylfield_controlblock *targetfield;
	IF1ARG(bgq_weylfield_controlblock *argfield1;)
	IF2ARG(bgq_weylfield_controlblock *argfield2;)
	IFEXTRA(NAME2(OPERATOR_NAME,extraparm_t) extra;)
} NAME2(OPERATOR_NAME,args_t);



static inline void NAME2(OPERATOR_NAME,worker)(void *arg_untyped, size_t tid, size_t threads, bool writeSloppy IF1ARG(, bool readWeyllayout1, bool sloppy1, bool mul1) IF2ARG(, bool readWeyllayout2, bool sloppy2, bool mul2)) {
	NAME2(OPERATOR_NAME,args_t) *arg = arg_untyped;
	tristate isOdd = arg->isOdd;
	bgq_weylfield_controlblock *targetfield = arg->targetfield;
	IF1ARG(bgq_weylfield_controlblock *argfield1 = arg->argfield1;)
	IF2ARG(bgq_weylfield_controlblock *argfield2 = arg->argfield2;)
	PREPROCESSOR_FOREACH( ,;,;, ,IDENTITY, OPERATOR_EXTRAPARMS)
	PREPROCESSOR_FOREACH(,,,,OPERATOR_EXTRA_ASSIGN, OPERATOR_EXTRAARGS)

	size_t workload = PHYSICAL_VOLUME;
	size_t threadload = (workload+threads-1)/threads;
	ucoord beginj = tid*threadload;
	ucoord endj = min_sizet((tid+1)*threadload,workload);

	IF1ARG(bgq_spinorfield_streamSpinor(argfield1, isOdd, beginj, readWeyllayout1, sloppy1, mul1, false);)
	IF2ARG(bgq_spinorfield_streamSpinor(argfield2, isOdd, beginj, readWeyllayout2, sloppy2, mul2, false);)
	for (ucoord ic = beginj; ic < endj; ic+=1) {
#ifndef NDEBUG
		assert(isOdd!=tri_unknown);
		ucoord ih = bgq_collapsed2halfvolume(isOdd, ic);
		ucoord tv = bgq_halfvolume2tv(ih);
		ucoord t1 = bgq_halfvolume2t1(isOdd, ih);
		ucoord t2 = bgq_halfvolume2t2(isOdd, ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);
#endif

#if OPERATOR_ARGFIELDS>=1
		bgq_su3_spinor_decl(spinor1);
		bgq_spinorfield_readSpinor(&spinor1, argfield1, isOdd, ic, readWeyllayout1, sloppy1, mul1, false);
		bgq_spinorfield_prefetchNextSpinor(argfield1, isOdd, ic, readWeyllayout1, sloppy1, mul1, false);
#endif

#if OPERATOR_ARGFIELDS>=2
		bgq_su3_spinor_decl(spinor2);
		bgq_spinorfield_readSpinor(&spinor2, argfield2, isOdd, ic, readWeyllayout2, sloppy2, mul2, false);
		bgq_spinorfield_prefetchNextSpinor(argfield2, isOdd, ic, readWeyllayout2, sloppy2, mul2, false);
#endif

		bgq_su3_spinor_decl(targetspinor);
		OPERATOR_VECSITEFUNC(bgq_su3_spinor_vars(&targetspinor) IF1ARG(, bgq_su3_spinor_vars(spinor1)) IF2ARG(, bgq_su3_spinor_vars(spinor2)) OPERATOR_EXTRAARGLIST, ic);

		// Warning: targetsite==argfield1 is possible and the inline assembler does not have memory barriers! If there is not data dependency, the compiler might arrange the stores before the loads!
		//REORDER_BARRIER
		if (writeSloppy) {
			bgq_spinorsite_float *targetsite = &targetfield->sec_fullspinor_float[ic];
			bgq_su3_spinor_store_float(targetsite, targetspinor);
		} else {
			bgq_spinorsite_double *targetsite = &targetfield->sec_fullspinor_double[ic];
			bgq_su3_spinor_store_double(targetsite, targetspinor);
		}
	}
}


#if OPERATOR_ARGFIELDS==0
#define OPERATOR_WORKERNAME(precision) NAME3(OPERATOR_NAME,worker,precision)

static void OPERATOR_WORKERNAME(double)(void *arg_untyped, size_t tid, size_t threads) {
	NAME2(OPERATOR_NAME,worker)(arg_untyped, tid, threads, false);
}

static void OPERATOR_WORKERNAME(float)(void *arg_untyped, size_t tid, size_t threads) {
	NAME2(OPERATOR_NAME,worker)(arg_untyped, tid, threads, true);
}

static bgq_worker_func NAME2(OPERATOR_NAME,worker_funcs)[2] =
	{ &OPERATOR_WORKERNAME(double), &OPERATOR_WORKERNAME(float) };


#elif OPERATOR_ARGFIELDS==1
#define OPERATOR_WORKERNAME(precision,layout1) NAME4(OPERATOR_NAME,worker,NAME2(OPERATOR_LAYOUTNAME,layout1),precision)
#define OPERATOR_WORKERINST(precision,layout1) \
		static void OPERATOR_WORKERNAME(precision,layout1)(void *arg_untyped, size_t tid, size_t threads) { \
			NAME2(OPERATOR_NAME,worker)(arg_untyped, tid, threads, NAME2(OPERATOR_PRECISION_ISSLOPPY,precision), (layout1)&ly_weyl, (layout1)&ly_sloppy, (layout1)&ly_mul); \
		}

OPERATOR_WORKERINST(double,ly_full_double)
OPERATOR_WORKERINST(double,ly_weyl_double)

OPERATOR_WORKERINST(float,ly_full_float)
OPERATOR_WORKERINST(float,ly_weyl_float)


static bgq_worker_func NAME2(OPERATOR_NAME,worker_funcs)[2][BGQ_SPINORFIELD_LAYOUT_COUNT] = {
	               /* 0=ly_full_double */           /* 1=ly_weyl_double */           /* 2=ly_full_float */           /* 3=ly_weyl_float */
	/* double */ { &OPERATOR_WORKERNAME(double,0), &OPERATOR_WORKERNAME(double,1), NULL,                           NULL },
	/* float  */ { NULL,                            NULL,                            &OPERATOR_WORKERNAME(float,2), &OPERATOR_WORKERNAME(float,3) }
};


#elif OPERATOR_ARGFIELDS==2
#define OPERATOR_WORKERNAME(precision,layout1,layout2) NAME5(OPERATOR_NAME,worker,NAME2(OPERATOR_LAYOUTNAME,layout1),NAME2(OPERATOR_LAYOUTNAME,layout2),precision)
#define OPERATOR_WORKERINST(precision,layout1,layout2) \
		static void OPERATOR_WORKERNAME(precision,layout1,layout2)(void *arg_untyped, size_t tid, size_t threads) { \
			NAME2(OPERATOR_NAME,worker)(arg_untyped, tid, threads, NAME2(OPERATOR_PRECISION_ISSLOPPY,precision), (layout1)&ly_weyl, (layout1)&ly_sloppy, (layout1)&ly_mul, (layout2)&ly_weyl, (layout2)&ly_sloppy, (layout2)&ly_mul); \
		}

OPERATOR_WORKERINST(double,ly_full_double,ly_full_double)
OPERATOR_WORKERINST(double,ly_weyl_double,ly_full_double)
OPERATOR_WORKERINST(double,ly_full_double,ly_weyl_double)
OPERATOR_WORKERINST(double,ly_weyl_double,ly_weyl_double)

OPERATOR_WORKERINST(float,ly_full_float,ly_full_float)
OPERATOR_WORKERINST(float,ly_weyl_float,ly_full_float)
OPERATOR_WORKERINST(float,ly_full_float,ly_weyl_float)
OPERATOR_WORKERINST(float,ly_weyl_float,ly_weyl_float)

static bgq_worker_func NAME2(OPERATOR_NAME,worker_funcs)[2][BGQ_SPINORFIELD_LAYOUT_COUNT][BGQ_SPINORFIELD_LAYOUT_COUNT] = {
/* double */
{
						 /* 0=ly_full_double */           /* 1=ly_weyl_double */           /* 2=ly_full_float */           /* 3=ly_weyl_float */
/* 0=ly_full_double */ { &OPERATOR_WORKERNAME(double,0,0), &OPERATOR_WORKERNAME(double,0,1), NULL, NULL },
/* 1=ly_weyl_double */ { &OPERATOR_WORKERNAME(double,1,0), &OPERATOR_WORKERNAME(double,1,1), NULL, NULL },
/* 2=ly_full_float */  { NULL,                              NULL,                              NULL, NULL },
/* 3=ly_weyl_float */  { NULL,                              NULL,                              NULL, NULL }
},
/* float */
{
						/* 0=ly_full_double */ /* 1=ly_weyl_double */  /* 2=ly_full_float */           /* 3=ly_weyl_float */
/* 0=ly_full_double */ { NULL, NULL, NULL,                             NULL },
/* 1=ly_weyl_double */ { NULL, NULL, NULL,                             NULL },
/* 2=ly_full_float */  { NULL, NULL, &OPERATOR_WORKERNAME(float,2,2), &OPERATOR_WORKERNAME(float,2,3) },
/* 3=ly_weyl_float */  { NULL, NULL, &OPERATOR_WORKERNAME(float,3,2), &OPERATOR_WORKERNAME(float,3,3) }
}};
#else
#error Number of field arguments not supported
#endif


static void NAME2(OPERATOR_NAME,selector)(bool sloppy, bgq_weylfield_controlblock *targetfield, tristate isOdd IF1ARG(, bgq_weylfield_controlblock *argfield1) IF2ARG(, bgq_weylfield_controlblock *argfield2) OPERATOR_EXTRAPARMLIST) {
	IF1ARG(isOdd = tristate_combine(isOdd, argfield1->isOdd);)
	IF2ARG(isOdd = tristate_combine(isOdd, argfield2->isOdd);)

#if OPERATOR_ARGFIELDS==0
	bgq_worker_func workerfunc = NAME2(OPERATOR_NAME,worker_funcs)[sloppy];
#elif OPERATOR_ARGFIELDS==1
	bgq_spinorfield_layout layout1 = bgq_spinorfield_prepareRead(argfield1, isOdd, true, true, true, true, true);
	if (!NAME2(OPERATOR_NAME,worker_funcs)[sloppy][layout1]) {
		layout1 = bgq_spinorfield_prepareRead(argfield1, isOdd, false, !sloppy, sloppy, false, false); // Force a rewrite
	}
	bgq_worker_func workerfunc = NAME2(OPERATOR_NAME,worker_funcs)[sloppy][layout1];
#elif OPERATOR_ARGFIELDS==2
	bgq_spinorfield_layout layout1 = bgq_spinorfield_prepareRead(argfield1, isOdd, true, true, true, true, true);
	bgq_spinorfield_layout layout2 = bgq_spinorfield_prepareRead(argfield2, isOdd, true, true, true, true, true);
	if (!NAME2(OPERATOR_NAME,worker_funcs)[sloppy][layout1][layout2]) {
		if (NAME2(OPERATOR_NAME,worker_funcs)[sloppy][(sloppy ? ly_full_float : ly_full_double)][layout2]) {
			layout1 = bgq_spinorfield_prepareRead(argfield1, isOdd, false, !sloppy, sloppy, false, false);
		} else if (NAME2(OPERATOR_NAME,worker_funcs)[sloppy][layout1][(sloppy ? ly_full_float : ly_full_double)]) {
			layout2 = bgq_spinorfield_prepareRead(argfield2, isOdd, false, !sloppy, sloppy, false, false);
		} else {
			// Force rewrite of both
			layout1 = bgq_spinorfield_prepareRead(argfield1, isOdd, false, !sloppy, sloppy, false, false);
			layout2 = bgq_spinorfield_prepareRead(argfield2, isOdd, false, !sloppy, sloppy, false, false);
		}
	}
	bgq_worker_func workerfunc = NAME2(OPERATOR_NAME,worker_funcs)[sloppy][layout1][layout2];
#endif
	assert(workerfunc);

	bgq_spinorfield_prepareWrite(targetfield, isOdd, sloppy ? ly_full_float : ly_full_double, false IF1ARG(|| targetfield==argfield1) IF2ARG(|| targetfield==argfield2));

	bgq_master_sync();
	static NAME2(OPERATOR_NAME,args_t) operator_args;
	NAME2(OPERATOR_NAME,args_t) call_args = {
			.isOdd = isOdd,
			.targetfield = targetfield
			IF1ARG(, .argfield1 = argfield1)
			IF2ARG(, .argfield2 = argfield2)
			IFEXTRA(, .extra = { OPERATOR_EXTRAARGS })
	};
	operator_args = call_args;
	bgq_master_call(workerfunc, &operator_args);
}


static void NAME2(OPERATOR_NAME,double)(bgq_weylfield_controlblock *targetfield, tristate isOdd IF1ARG(, bgq_weylfield_controlblock *argfield1) IF2ARG(, bgq_weylfield_controlblock *argfield2) OPERATOR_EXTRAPARMLIST) {
	NAME2(OPERATOR_NAME,selector)(false, targetfield, isOdd IF1ARG(, argfield1) IF2ARG(, argfield2) OPERATOR_EXTRAARGLIST);
}


static void NAME2(OPERATOR_NAME,float)(bgq_weylfield_controlblock *targetfield, tristate isOdd IF1ARG(, bgq_weylfield_controlblock *argfield1) IF2ARG(, bgq_weylfield_controlblock *argfield2) OPERATOR_EXTRAPARMLIST) {
	NAME2(OPERATOR_NAME,selector)(true, targetfield, isOdd IF1ARG(, argfield1) IF2ARG(, argfield2) OPERATOR_EXTRAARGLIST);
}

#undef OPERATOR_EXTRAARGLIST
#undef OPERATOR_EXTRAPARMLIST
#undef OPERATOR_WORKERNAME
#undef OPERATOR_WORKERINST

#undef IFNOARG
#undef IF1ARG
#undef IF2ARG
#undef IFEXTRA

#undef OPERATOR_NAME
#undef OPERATOR_OUTPLACENAME
#undef OPERATOR_INPLACENAME
#undef OPERATOR_ARGFIELDS
#undef OPERATOR_EXTRAPARMS
#undef OPERATOR_EXTRAARGS
#undef OPERATOR_VECSITEFUNC




void bgq_spinorfield_rmul_plain_sub_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield1, bgq_weylfield_controlblock *sourcefield2, double r) {
	bgq_vector4double_decl(qr);
	bgq_cconst(qr, r, r);
	bgq_spinorfield_rmul_plain_sub_raw_double(targetfield, isOdd, sourcefield1, sourcefield2, bgq_vars(qr));
}


