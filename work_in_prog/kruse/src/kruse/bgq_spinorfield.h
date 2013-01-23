/*
 * bgq_spinorfield.h
 *
 *  Created on: Oct 26, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_SPINORFIELD_H_
#define BGQ_SPINORFIELD_H_

#include "bgq_field.h"
#include "bgq_utils.h"
#include "bgq_qpx.h"

#include <stdbool.h>

#ifndef BGQ_SPINORFIELD_C_
#define EXTERN_INLINE EXTERN_INLINE_DECLARATION
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE EXTERN_INLINE_DEFINITION
#define EXTERN_FIELD
#define EXTERN_INIT(val) = (val)
#endif


typedef struct {
	COMPLEX_PRECISION c[3];
} bgq_su3vector;

typedef struct {
	bgq_su3vector v[4];
} bgq_spinor;
typedef bgq_spinor bgq_spinor_nonvec;

typedef struct {
	COMPLEX_PRECISION s[2][3]; // 96 byte
} bgq_weyl_nonvec;


typedef enum {
	ly_weyl = (1<<0),
	ly_sloppy = (1<<1),
	ly_mul = (1<<2)
} bgq_spinorfield_layoutflags;
typedef enum {
	ly_full_double=0,
	ly_weyl_double=ly_weyl,
	ly_full_float=ly_sloppy,
	ly_weyl_float=ly_weyl|ly_sloppy,
	ly_full_double_mul=ly_mul,
	ly_weyl_double_mul=ly_mul|ly_weyl,
	ly_full_float_mul=ly_mul|ly_sloppy,
	ly_weyl_float_mul=ly_mul|ly_weyl|ly_sloppy,
	ly_legacy=(1<<3),
	ly_none=-1,
} bgq_spinorfield_layout;
#define BGQ_SPINORFIELD_LAYOUT_COUNT 9




typedef struct {
	//bool isInitialized;
	tristate isOdd;
	//bool ptr_isOdd;
	//bool hasWeylfieldData;
	//bool isWeyllayoutSloppy;
	//bool waitingForRecv; /* true==Need to wait for SPI recv and then copy data to consecutive area; false==All data available in sec_surface and sec_body */
	//bool waitingForRecvNoSPI;
	bgq_hmflags hmflags;
	bool pendingDatamove;
	//bool hasFullspinorData;
	//bool isFulllayoutSloppy;

	spinor *legacy_field;
	bool has_legacy;

	bool has_weyllayout_double;
	bool has_weyllayout_float;
	bool has_fulllayout_double;
	bool has_fulllayout_float;

	//uint8_t *sec_weyl;
	//bgq_weyl_vec *sec_index; // obsolete
	//bgq_weyl_vec *sec_send[PHYSICAL_LD]; // obsolete
	//bgq_weyl_vec *sec_recv[PHYSICAL_LD]; // obsolete
	bgq_weylsite_double *sec_collapsed_double;
	bgq_weylsite_float *sec_collapsed_float;
	//bgq_weylsite *sec_surface;
	//bgq_weylsite *sec_body;
	//uint8_t *sec_end;


	bgq_spinorsite_double *sec_fullspinor_double;
	//bgq_spinorsite *sec_fullspinor_surface;
	//bgq_spinorsite *sec_fullspinor_body;
	bgq_spinorsite_float *sec_fullspinor_float;


	//TODO: We may even interleave these with the data itself, but may cause alignment issues
	// Idea: sizeof(bgq_weyl_ptr_t)==10*8==80, so one bgq_weyl_ptr_t every 2(5;10) spinors solves the issue
	// In case we write as fullspinor layout, the are not needed
	bgq_weyl_ptr_t_double *sendptr_double[PHYSICAL_LP];
	bgq_weyl_vec_double **consptr_double[PHYSICAL_LP][PHYSICAL_LD];

	bgq_weyl_ptr_t_float *sendptr_float[PHYSICAL_LP];
	bgq_weyl_vec_float **consptr_float[PHYSICAL_LP][PHYSICAL_LD];

	struct bgq_weylfield_collection *collectionBase;
} bgq_weylfield_controlblock;


#define BGQ_SEC_FULLLAYOUT NAME2(sec_fullspinor,PRECISION)
#define BGQ_SEC_WEYLLAYOUT NAME2(sec_collapsed,PRECISION)
#define BGQ_HAS_FULLLAYOUT NAME2(has_fulllayout,PRECISION)
#define BGQ_HAS_WEYLLAYOUT NAME2(has_collapsed,PRECISION)
#define BGQ_SENDPTR NAME2(sendptr,PRECISION)
#define BGQ_CONSPTR NAME2(consptr,PRECISION)


EXTERN_FIELD bgq_weylfield_controlblock *g_bgq_spinorfields EXTERN_INIT(NULL);

typedef struct bgq_weylfield_collection {
	const spinor *legacy_base;
	size_t fieldsize;
	size_t count;
	struct bgq_weylfield_collection *prev;
	struct bgq_weylfield_collection *next;
	bgq_weylfield_controlblock controlblocks[];
} bgq_weylfield_collection;


void bgq_spinorfields_init(void);
bgq_weylfield_collection *bgq_spinorfields_allocate(size_t count, spinor *legacyFields, size_t fieldLength);
void bgq_spinorfields_free(bgq_weylfield_collection *collection);

size_t bgq_pointer2offset_raw(bgq_weylfield_controlblock *field, void *ptr, bool check);
size_t bgq_pointer2offset(bgq_weylfield_controlblock *field, void *ptr);




#define bgq_weyl_fromqpx(arg) bgq_weyl_fromqpx_raw(bgq_su3_weyl_vars(arg))
EXTERN_INLINE bgq_weyl_vec_double bgq_weyl_fromqpx_raw(bgq_su3_weyl_params(weyl)) {
	bgq_weyl_vec_double result;
	result.s[0][0][0] = bgq_cmplxval1(weyl_v0_c0);
	result.s[0][0][1] = bgq_cmplxval2(weyl_v0_c0);
	result.s[0][1][0] = bgq_cmplxval1(weyl_v0_c1);
	result.s[0][1][1] = bgq_cmplxval2(weyl_v0_c1);
	result.s[0][2][0] = bgq_cmplxval1(weyl_v0_c2);
	result.s[0][2][1] = bgq_cmplxval2(weyl_v0_c2);
	result.s[1][0][0] = bgq_cmplxval1(weyl_v1_c0);
	result.s[1][0][1] = bgq_cmplxval2(weyl_v1_c0);
	result.s[1][1][0] = bgq_cmplxval1(weyl_v1_c1);
	result.s[1][1][1] = bgq_cmplxval2(weyl_v1_c1);
	result.s[1][2][0] = bgq_cmplxval1(weyl_v1_c2);
	result.s[1][2][1] = bgq_cmplxval2(weyl_v1_c2);
	return result;
}


#define bgq_weyl_fromqpxk(arg,k) bgq_weyl_fromqpxk_raw(bgq_su3_weyl_vars(arg), k)
EXTERN_INLINE bgq_weyl_nonvec bgq_weyl_fromqpxk_raw(bgq_su3_weyl_params(weyl), ucoord k) {
	bgq_weyl_nonvec result;
	result.s[0][0] = bgq_cmplxval(weyl_v0_c0,k);
	result.s[0][1] = bgq_cmplxval(weyl_v0_c1,k);
	result.s[0][2] = bgq_cmplxval(weyl_v0_c2,k);
	result.s[1][0] = bgq_cmplxval(weyl_v1_c0,k);
	result.s[1][1] = bgq_cmplxval(weyl_v1_c1,k);
	result.s[1][2] = bgq_cmplxval(weyl_v1_c2,k);
	return result;
}



EXTERN_INLINE bgq_weyl_nonvec bgq_weyl_extractvec(bgq_weyl_vec_double weylvec, size_t k) {
	assert(0 <= k && k < PHYSICAL_LK);
	bgq_weyl_nonvec result;
	for (size_t v = 0; v < 2; v+=1) {
		for (size_t c = 0; c < 3; c+=1) {
			result.s[v][c] = weylvec.s[v][c][k];
		}
	}
	return result;
}


#define bgq_spinor_fromqpx(spinor,k) bgq_spinor_fromqpx_raw(bgq_su3_spinor_vars(spinor),k)
EXTERN_INLINE bgq_spinor bgq_spinor_fromqpx_raw(bgq_su3_spinor_params(spinor), ucoord k) {
	assert(0 <= k && k < PHYSICAL_LK);

	bgq_spinor result;
	result.v[0].c[0] = bgq_cmplxval(spinor_v0_c0,k);
	result.v[0].c[1] = bgq_cmplxval(spinor_v0_c1,k);
	result.v[0].c[2] = bgq_cmplxval(spinor_v0_c2,k);
	result.v[1].c[0] = bgq_cmplxval(spinor_v1_c0,k);
	result.v[1].c[1] = bgq_cmplxval(spinor_v1_c1,k);
	result.v[1].c[2] = bgq_cmplxval(spinor_v1_c2,k);
	result.v[2].c[0] = bgq_cmplxval(spinor_v2_c0,k);
	result.v[2].c[1] = bgq_cmplxval(spinor_v2_c1,k);
	result.v[2].c[2] = bgq_cmplxval(spinor_v2_c2,k);
	result.v[3].c[0] = bgq_cmplxval(spinor_v3_c0,k);
	result.v[3].c[1] = bgq_cmplxval(spinor_v3_c1,k);
	result.v[3].c[2] = bgq_cmplxval(spinor_v3_c2,k);
	return result;
}


EXTERN_INLINE bgq_spinor bgq_spinor_fromlegacy(spinor *sp) {
	bgq_spinor result;
	result.v[0].c[0] = sp->s0.c0;
	result.v[0].c[1] = sp->s0.c1;
	result.v[0].c[2] = sp->s0.c2;
	result.v[1].c[0] = sp->s1.c0;
	result.v[1].c[1] = sp->s1.c1;
	result.v[1].c[2] = sp->s1.c2;
	result.v[2].c[0] = sp->s2.c0;
	result.v[2].c[1] = sp->s2.c1;
	result.v[2].c[2] = sp->s2.c2;
	result.v[3].c[0] = sp->s3.c0;
	result.v[3].c[1] = sp->s3.c1;
	result.v[3].c[2] = sp->s3.c2;
	return result;
}

#define bgq_spinor_fromvec NAME2(bgq_spinor_fromvec,PRECISION)
EXTERN_INLINE bgq_spinor bgq_spinor_fromvec_double(bgq_spinor_vec_double spinorvec, ucoord k) {
	assert(0 <= k && k < PHYSICAL_LK);

	bgq_spinor result;
	for (ucoord v = 0; v < 4; v+=1) {
		for (ucoord c = 0; c < 3; c+=1) {
			result.v[v].c[c] = spinorvec.s[v][c][k];
		}
	}
	return result;
}
EXTERN_INLINE bgq_spinor bgq_spinor_fromvec_float(bgq_spinor_vec_float spinorvec, ucoord k) {
	assert(0 <= k && k < PHYSICAL_LK);

	bgq_spinor result;
	for (ucoord v = 0; v < 4; v+=1) {
		for (ucoord c = 0; c < 3; c+=1) {
			result.v[v].c[c] = spinorvec.s[v][c][k];
		}
	}
	return result;
}


#define bgq_weyl_extractfromqpxvec(weyl, k) bgq_weyl_extractfromqpxvec_raw(bgq_su3_weyl_vars(weyl),k)
EXTERN_INLINE bgq_weyl_nonvec bgq_weyl_extractfromqpxvec_raw(bgq_su3_weyl_params(weyl), size_t k) {
	bgq_weyl_vec_double vec = bgq_weyl_fromqpx(weyl);
	return bgq_weyl_extractvec(vec,k);
}



//void bgq_spinorfield_enableLayout(bgq_weylfield_controlblock *field, bool isOdd, bgq_spinorfield_layout layout, bool disableOthers);
//void bgq_spinorfield_setup(bgq_weylfield_controlblock *field, bool isOdd, bool readFullspinor, bool writeFullspinor, bool readWeyl, bool writeWeyl, bool writeFloat); // obsolete
//void bgq_spinorfield_setup_double(bgq_weylfield_controlblock *field, bool isOdd, bool readFullspinor, bool writeFullspinor, bool readWeyl, bool writeWeyl); // obsolete
//void bgq_spinorfield_setup_float(bgq_weylfield_controlblock *field, bool isOdd, bool readFullspinor, bool writeFullspinor, bool readWeyl, bool writeWeyl); // obsolete
//void bgq_spinorfield_transfer(bool isOdd, bgq_weylfield_controlblock *targetfield, spinor* sourcefield);
double bgq_spinorfield_compare(bool isOdd, bgq_weylfield_controlblock *bgqfield, bgq_weylfield_controlblock *reffield, bool silent);

typedef struct {
	bgq_weylfield_controlblock *spinorfield;
	bgq_hmflags opts;
} bgq_work_datamove;
void bgq_HoppingMatrix_datamovet_worker(void *arg_untyped, size_t tid, size_t threads);



EXTERN_INLINE void bgq_spinor_expect(bgq_spinor spinor, scoord t, scoord x, scoord y, scoord z) {
#ifdef BGQ_COORDCHECK
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	ucoord t_global = bgq_local2global_t(t);
	ucoord x_global = bgq_local2global_x(x);
	ucoord y_global = bgq_local2global_y(y);
	ucoord z_global = bgq_local2global_z(z);

	assert(spinor.v[0].c[0] == t_global);
	assert(spinor.v[0].c[1] == x_global);
	assert(spinor.v[0].c[2] == y_global);
	assert(spinor.v[1].c[0] == z_global);
	assert(spinor.v[1].c[1] == 0);
	assert(spinor.v[1].c[2] == 0);
	assert(spinor.v[2].c[0] == 0);
	assert(spinor.v[2].c[1] == 0);
	assert(spinor.v[2].c[2] == 0);
	assert(spinor.v[3].c[0] == 0);
	assert(spinor.v[3].c[1] == 0);
	assert(spinor.v[3].c[2] == 0);
#endif
}


#define bgq_spinorveck_expect NAME2(bgq_spinorveck_expect,PRECISION)

#ifdef BGQ_COORDCHECK
#define bgq_spinorveck_expect_double(spinor,k,t,x,y,z) bgq_spinorveck_expect_double_raw(bgq_su3_spinor_vars(spinor),k,t,x,y,z)
#else
#define bgq_spinorveck_expect_double(spinor,k,t,x,y,z)
#endif
EXTERN_INLINE void bgq_spinorveck_expect_double_raw(bgq_spinor_vec_double spinor, ucoord k, scoord t,scoord x,scoord y,scoord z) {
	assert(0 <= k && k < PHYSICAL_LK);
	bgq_spinor_expect(bgq_spinor_fromvec_double(spinor,k),t,x,y,z);
}


#ifdef BGQ_COORDCHECK
#define bgq_spinorveck_expect_float(spinor,k,t,x,y,z) bgq_spinorveck_expect_float_raw(bgq_su3_spinor_vars(spinor),k,t,x,y,z)
#else
#define bgq_spinorveck_expect_float(spinor,k,t,x,y,z)
#endif
EXTERN_INLINE void bgq_spinorveck_expect_float_raw(bgq_spinor_vec_float spinor, ucoord k, scoord t,scoord x,scoord y,scoord z) {
	assert(0 <= k && k < PHYSICAL_LK);
	bgq_spinor_expect(bgq_spinor_fromvec_float(spinor,k),t,x,y,z);
}


EXTERN_INLINE void bgq_spinorvec_expect(bgq_spinor_vec_double spinor, scoord t1, scoord t2, scoord x,scoord y,scoord z) {
	bgq_spinor_expect(bgq_spinor_fromvec_double(spinor,0),t1,x,y,z);
	bgq_spinor_expect(bgq_spinor_fromvec_double(spinor,1),t2,x,y,z);
}


#ifdef BGQ_COORDCHECK
#define bgq_spinorqpx_expect(spinor,t_left,t_right,x,y,z) bgq_spinorqpx_expect_raw(bgq_su3_spinor_vars(spinor),t_left,t_right,x,y,z)
#else
#define bgq_spinorqpx_expect(spinor,t_left,t_right,x,y,z)
#endif
EXTERN_INLINE void bgq_spinorqpx_expect_raw(bgq_su3_spinor_params(spinor),scoord t_left,scoord t_right,scoord x,scoord y,scoord z) {
	bgq_spinor_expect(bgq_spinor_fromqpx(spinor,0),t_left,x,y,z);
	bgq_spinor_expect(bgq_spinor_fromqpx(spinor,1),t_right,x,y,z);
}


#ifdef BGQ_COORDCHECK
#define bgq_spinorlegacy_expect(spinor,t,x,y,z) bgq_spinorlegacy_expect_raw(spinor,t,x,y,z)
#else
#define bgq_spinorlegacy_expect(spinor,t,x,y,z)
#endif
EXTERN_INLINE void bgq_spinorlegacy_expect_raw(spinor *sp,scoord t,scoord x,scoord y,scoord z) {
	bgq_spinor_expect(bgq_spinor_fromlegacy(sp),t,x,y,z);
}

#ifdef BGQ_COORDCHECK
#define bgq_spinorqpxk_expect(spinor,k,t,x,y,z) bgq_spinorqpxk_expect_raw(bgq_su3_spinor_vars(spinor),k,t,x,y,z)
#else
#define bgq_spinorqpxk_expect(spinor,k,t,x,y,z)
#endif
EXTERN_INLINE void bgq_spinorqpxk_expect_raw(bgq_su3_spinor_params(spinor),ucoord k,scoord t,scoord x,scoord y,scoord z) {
	bgq_spinor_expect(bgq_spinor_fromqpx(spinor,k),t,x,y,z);
}


void bgq_weyl_expect(bgq_weyl_nonvec weyl, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc);

#ifdef BGQ_COORDCHECK
#define bgq_weylqpx_expect(weyl,t1,t2,x,y,z,d,isSrc) bgq_weylqpx_expect_raw(bgq_su3_weyl_vars(weyl),t1,t2,x,y,z,d,isSrc)
#else
#define bgq_weylqpx_expect(weyl,t1,t2,x,y,z,d,isSrc)
#endif
EXTERN_INLINE void bgq_weylqpx_expect_raw(bgq_su3_weyl_params(weyl), ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
	bgq_weyl_expect(bgq_weyl_fromqpxk(weyl,0),t1,x,y,z,d,isSrc);
	bgq_weyl_expect(bgq_weyl_fromqpxk(weyl,1),t2,x,y,z,d,isSrc);
}


#ifdef BGQ_COORDCHECK
#define bgq_weylqpxk_expect(weyl,k,t,x,y,z,d,isSrc) bgq_weylqpxk_expect_raw(bgq_su3_weyl_vars(weyl),k,t,x,y,z,d,isSrc)
#else
#define bgq_weylqpxk_expect(weyl, k, t, x, y, z, d, isSrc)
#endif
EXTERN_INLINE void bgq_weylqpxk_expect_raw(bgq_su3_weyl_params(weyl), ucoord k, ucoord t,  ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
	bgq_weyl_expect(bgq_weyl_fromqpxk(weyl, k), t, x, y, z, d, isSrc);
}


#ifdef BGQ_COORDCHECK
#define BGQ_WEYLLEFT_EXPECT(weyl,t,x,y,z) \
		bgq_spinorfield_weyl_expect(bgq_weyl_extractfromqpxvec(weyl,0),t,x,y,z)
#define BGQ_WEYLRIGHT_EXPECT(weyl,t,x,y,z) \
		bgq_spinorfield_weyl_expect(bgq_weyl_extractfromqpxvec(weyl,1),t,x,y,z)
#define BGQ_WEYL_EXPECT(weyl,t1,t2,x,y,z) \
		BGQ_WEYLLEFT_EXPECT(weyl,t1,x,y,z); \
		BGQ_WEYLRIGHT_EXPECT(weyl,t2,x,y,z)
#else
#define BGQ_WEYLLEFT_EXPECT(weyl,t,x,y,z)
#define BGQ_WEYLRIGHT_EXPECT(weyl,t,x,y,z)
#define BGQ_WEYL_EXPECT(weyl,t,x,y,z)
#endif


EXTERN_INLINE bgq_weyl_nonvec bgq_weyl_fromvec(bgq_weyl_vec_double weylvec, ucoord k) {
	assert(0 <= k && k < PHYSICAL_LK);
	bgq_weyl_nonvec result;
	for (size_t v = 0; v < 2; v += 1) {
		for (size_t c = 0; c < 3; c += 1) {
			result.s[v][c] = weylvec.s[v][c][k];
		}
	}
	return result;
}


void bgq_weyl_expect(bgq_weyl_nonvec weyl, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc);

EXTERN_INLINE void bgq_weylvec_expect(bgq_weyl_vec_double weyl, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
#ifdef BGQ_COORDCHECK
	bgq_weyl_expect(bgq_weyl_fromvec(weyl, 0), t1, x, y, z, d, isSrc);
	bgq_weyl_expect(bgq_weyl_fromvec(weyl, 1), t2, x, y, z, d, isSrc);
#endif
}



bgq_spinor bgq_legacy_getspinor(spinor *spinor, ucoord t, ucoord x, ucoord y, ucoord z);
bgq_spinor bgq_spinorfield_getspinor(bgq_weylfield_controlblock *field, size_t t, size_t x, size_t y, size_t z) ;

void bgq_weylveck_written_double(bgq_weyl_vec_double *targetweyl, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc);
void bgq_weylveck_written_float(bgq_weyl_vec_float *targetweyl, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc);


#define bgq_weylvec_written NAME2(bgq_weylvec_written,PRECISION)

#ifdef BGQ_COORDCHECK
#define bgq_weylvec_written_double(targetweyl, t1, t2, x, y, z, d, isSrc) bgq_weylvec_written_impl_double(targetweyl, t1, t2, x, y, z, d, isSrc)
#else
#define bgq_weylvec_written_double(targetweyl, t1, t2, x, y, z, d, isSrc)
#endif
EXTERN_INLINE void bgq_weylvec_written_impl_double(bgq_weyl_vec_double *targetweyl, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
	bgq_weylveck_written_double(targetweyl, 0, t1, x, y, z, d, isSrc);
	bgq_weylveck_written_double(targetweyl, 1, t2, x, y, z, d, isSrc);
}

#ifdef BGQ_COORDCHECK
#define bgq_weylvec_written_float(targetweyl, t1, t2, x, y, z, d, isSrc) bgq_weylvec_written_impl_float(targetweyl, t1, t2, x, y, z, d, isSrc)
#else
#define bgq_weylvec_written_float(targetweyl, t1, t2, x, y, z, d, isSrc)
#endif
EXTERN_INLINE void bgq_weylvec_written_impl_float(bgq_weyl_vec_float *targetweyl, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
	bgq_weylveck_written_float(targetweyl, 0, t1, x, y, z, d, isSrc);
	bgq_weylveck_written_float(targetweyl, 1, t2, x, y, z, d, isSrc);
}

void bgq_spinorveck_written_double(bgq_spinorsite_double *targetspinor, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z);
void bgq_spinorveck_written_float(bgq_spinorsite_float *targetspinor, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z);

#ifdef BGQ_COORDCHECK
#define bgq_spinorvec_written_double(target, t1, t2, x, y, z) bgq_spinorvec_written_impl_double(target, t1, t2, x, y, z)
#else
#define bgq_spinorvec_written_double(target, t1, t2, x, y, z)
#endif
EXTERN_INLINE void bgq_spinorvec_written_impl_double(bgq_spinor_vec_double *target, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z) {
	bgq_spinorveck_written_double(target, 0, t1, x, y, z);
	bgq_spinorveck_written_double(target, 1, t2, x, y, z);
}

#ifdef BGQ_COORDCHECK
#define bgq_spinorvec_written_float(target, t1, t2, x, y, z) bgq_spinorvec_written_impl_float(target, t1, t2, x, y, z)
#else
#define bgq_spinorvec_written_float(target, t1, t2, x, y, z)
#endif
EXTERN_INLINE void bgq_spinorvec_written_impl_float(bgq_spinor_vec_float *target, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z) {
	bgq_spinorveck_written_float(target, 0, t1, x, y, z);
	bgq_spinorveck_written_float(target, 1, t2, x, y, z);
}



typedef enum {
	BGQREF_SOURCE,

	BGQREF_TDOWN_SOURCE,
	BGQREF_TDOWN,
	BGQREF_TDOWN_GAUGE,
	BGQREF_TDOWN_WEYL,
	BGQREF_TDOWN_KAMUL,

	BGQREF_TUP_SOURCE,
	BGQREF_TUP,
	BGQREF_TUP_GAUGE,
	BGQREF_TUP_WEYL,
	BGQREF_TUP_KAMUL,

	BGQREF_XUP,
	BGQREF_XUP_GAUGE,
	BGQREF_XUP_WEYL,
	BGQREF_XUP_KAMUL,

	BGQREF_XDOWN,
	BGQREF_XDOWN_GAUGE,
	BGQREF_XDOWN_WEYL,
	BGQREF_XDOWN_KAMUL,

	BGQREF_ZUP,
	BGQREF_ZUP_GAUGE,
	BGQREF_ZUP_KAMUL,


	BGQREF_TUP_RECV,

	BGQREF_TDOWN_RECV,
	BGQREF_TDOWN_ACCUM,

	BGQREF_XUP_RECV,
	BGQREF_XUP_ACCUM,

	BGQREF_XDOWN_RECV,
	BGQREF_XDOWN_ACCUM,

	BGQREF_YUP_RECV,
	BGQREF_YUP_ACCUM,

	BGQREF_YDOWN_ACCUM,

	BGQREF_ZUP_RECV,
	BGQREF_ZUP_ACCUM,

	BGQREF_ACCUM,
	BGQREF_RESULT,

	BGQREF_count
} bgqref;


#ifndef BGQ_REFCOMPARE
#define bgq_initbgqref()
#define bgq_setdesc(idx,desc)
#define bgq_setrefvalue( t,  x,  y,  z,  idx,  val)
#define bgq_setbgqvalue( t,  x,  y,  z,  idx,  val)
#define bgq_setbgqvalue_src( t,  x,  y,  z,  d,  idx,  val)
#define bgq_savebgqref()
#else
#define bgq_initbgqref() bgq_initbgqref_impl()
#define bgq_setdesc(idx,desc) bgq_setdesc_impl(idx,desc)
#define bgq_setrefvalue( t,  x,  y,  z,  idx,  val) bgq_setrefvalue_impl( t,  x,  y,  z,  idx,  val)
#define bgq_setbgqvalue( t,  x,  y,  z,  idx,  val) bgq_setbgqvalue_impl( t,  x,  y,  z,  idx,  val)
#define bgq_setbgqvalue_src( t,  x,  y,  z,  d,  idx,  val) bgq_setbgqvalue_src_impl( t,  x,  y,  z,  d,  idx,  val)
#define bgq_savebgqref() bgq_savebgqref_impl()
#endif

void bgq_initbgqref_impl(void);
void bgq_setdesc_impl(int idx, char *desc);
void bgq_setrefvalue_impl(int t, int x, int y, int z, bgqref idx, complexdouble val);
void bgq_setbgqvalue_impl(int t, int x, int y, int z, bgqref idx, complexdouble val);
void bgq_setbgqvalue_src_impl(ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bgqref idx, complexdouble val);
void bgq_savebgqref_impl(void);

//size_t bgq_fieldpointer2offset(void *ptr);
bgq_weyl_vec_double *bgq_section_baseptr_double(bgq_weylfield_controlblock *field, bgq_weylfield_section section);
bgq_weyl_vec_float *bgq_section_baseptr_float(bgq_weylfield_controlblock *field, bgq_weylfield_section section);



EXTERN_INLINE void bgq_spinorfield_streamSpinor(bgq_weylfield_controlblock *field, tristate isOdd, ucoord ic_begin, bool readWeyllayout, bool sloppy, bool mul, bool legacy) {
	if (legacy) {
		// No prefetch for legacy field
	} else if (readWeyllayout) {
		if (sloppy) {
			assert(field->sec_collapsed_float);
			bgq_prefetch_forward(&field->sec_collapsed_float[ic_begin]);
		} else {
			assert(field->sec_collapsed_double);
			bgq_prefetch_forward(&field->sec_collapsed_double[ic_begin]);
		}
	} else {
		if (sloppy) {
			assert(field->sec_fullspinor_float);
			bgq_prefetch_forward(&field->sec_fullspinor_float[ic_begin]);
		} else {
			assert(field->sec_fullspinor_double);
			bgq_prefetch_forward(&field->sec_fullspinor_double[ic_begin]);
		}
	}
}


EXTERN_INLINE void bgq_spinorfield_prefetchSpinor(bgq_weylfield_controlblock *field, tristate isOdd, ucoord ic, bool readWeyllayout, bool sloppy, bool mul, bool legacy) {
	if (legacy) {
		// No prefetch for legacy field
	} else if (readWeyllayout) {
		if (sloppy) {
			assert(field->sec_collapsed_float);
			// Prefetch only first element; the next ones should be prefetched by the reading function
			bgq_su3_weyl_prefetch_float(&field->sec_collapsed_float[ic].d[0]);
		} else {
			assert(field->sec_collapsed_double);
			// Prefetch only first element; the next ones should be prefetched by the reading function
			bgq_su3_weyl_prefetch_float(&field->sec_collapsed_double[ic].d[0]);
		}
	} else {
		if (sloppy) {
			assert(field->sec_fullspinor_float);
			bgq_su3_spinor_prefetch_float(&field->sec_fullspinor_float[ic]);
		} else {
			assert(field->sec_fullspinor_double);
			bgq_su3_spinor_prefetch_float(&field->sec_fullspinor_double[ic]);
		}
	}
}


EXTERN_INLINE void bgq_spinorfield_prefetchNextSpinor(bgq_weylfield_controlblock *field, tristate isOdd, ucoord ic, bool readWeyllayout, bool sloppy, bool mul, bool legacy) {
	if (legacy) {
		// No prefetch for legacy field
	} else if (readWeyllayout) {
		if (sloppy) {
			assert(field->sec_collapsed_float);
			// Prefetch only first element; the next ones should be prefetched by the reading function
			bgq_su3_weyl_prefetch_float(&field->sec_collapsed_float[ic+1].d[0]);
		} else {
			assert(field->sec_collapsed_double);
			// Prefetch only first element; the next ones should be prefetched by the reading function
			bgq_su3_weyl_prefetch_float(&field->sec_collapsed_double[ic+1].d[0]);
		}
	} else {
		if (sloppy) {
			assert(field->sec_fullspinor_float);
			bgq_su3_spinornext_prefetch_float(&field->sec_fullspinor_float[ic]);
		} else {
			assert(field->sec_fullspinor_double);
			bgq_su3_spinornext_prefetch_float(&field->sec_fullspinor_double[ic]);
		}
	}
}


#define bgq_spinorfield_readSpinor(target, field, isOdd, ic, readWeyllayout, sloppy, mul, legacy) bgq_spinorfield_readSpinor_raw(bgq_su3_spinor_vars(target), field, isOdd, ic, readWeyllayout, sloppy, mul, legacy)
EXTERN_INLINE void bgq_spinorfield_readSpinor_raw(bgq_su3_spinor_params(*target), bgq_weylfield_controlblock *field, tristate isOdd, ucoord ic, bool readWeyllayout, bool sloppy, bool mul, bool legacy) {
	assert(field);
	//assert(field->isInitialized);

#ifndef NDEBUG
	ucoord ih = bgq_collapsed2halfvolume(isOdd, ic);
	ucoord t1 = bgq_halfvolume2t1(isOdd, ih);
	ucoord t2 = bgq_halfvolume2t2(isOdd, ih);
	ucoord tv = bgq_halfvolume2tv(ih);
	ucoord x = bgq_halfvolume2x(ih);
	ucoord y = bgq_halfvolume2y(ih);
	ucoord z = bgq_halfvolume2z(ih);
#endif

	bgq_su3_spinor_decl(spinor);
	if (legacy) {
		assert(isOdd != tri_unknown);
		int i1_eosub = bgq_collapsed2eosub(isOdd, ic, 0);
		spinor *addr1 = &field->legacy_field[i1_eosub];
		bgq_vector4double_decl(left0);
		bgq_vector4double_decl(left1);
		bgq_vector4double_decl(left2);
		bgq_vector4double_decl(left3);
		bgq_vector4double_decl(left4);
		bgq_vector4double_decl(left5);
		bgq_qvlfduxa(left0, addr1, 0);
		bgq_qvlfduxa(left1, addr1, 32);
		bgq_qvlfduxa(left2, addr1, 32);
		bgq_qvlfduxa(left3, addr1, 32);
		bgq_qvlfduxa(left4, addr1, 32);
		bgq_qvlfduxa(left5, addr1, 32);

		int i2_eosub = bgq_collapsed2eosub(isOdd, ic, 1);
		spinor *addr2 = &field->legacy_field[i2_eosub];
		bgq_vector4double_decl(right0);
		bgq_vector4double_decl(right1);
		bgq_vector4double_decl(right2);
		bgq_vector4double_decl(right3);
		bgq_vector4double_decl(right4);
		bgq_vector4double_decl(right5);
		bgq_qvlfduxa(right0, addr2, 0);
		bgq_qvlfduxa(right1, addr2, 32);
		bgq_qvlfduxa(right2, addr2, 32);
		bgq_qvlfduxa(right3, addr2, 32);
		bgq_qvlfduxa(right4, addr2, 32);
		bgq_qvlfduxa(right5, addr2, 32);

		bgq_lmerge(spinor_v0_c0, left0, right0);
		bgq_rmerge(spinor_v0_c1, left0, right0);
		bgq_lmerge(spinor_v0_c2, left1, right1);
		bgq_rmerge(spinor_v1_c0, left1, right1);
		bgq_lmerge(spinor_v1_c1, left2, right2);
		bgq_rmerge(spinor_v1_c2, left2, right2);
		bgq_lmerge(spinor_v2_c0, left3, right3);
		bgq_rmerge(spinor_v2_c1, left3, right3);
		bgq_lmerge(spinor_v2_c2, left4, right4);
		bgq_rmerge(spinor_v3_c0, left4, right4);
		bgq_lmerge(spinor_v3_c1, left5, right5);
		bgq_rmerge(spinor_v3_c2, left5, right5);
	} else {
		if (readWeyllayout) {
			//assert(field->hasWeylfieldData);
			if (sloppy) {
				//assert(field->isWeyllayoutSloppy);
				bgq_weylsite_float *weylsite = (bgq_weylsite_float *)((uint8_t*)&field->sec_collapsed_float[ic] - 16);

				#define PRECISION float
				#define BGQ_READWEYLLAYOUT_INC_
				#include "bgq_ReadWeyllayout.inc.c"
				#undef PRECISION
			} else {
				//assert(!field->isWeyllayoutSloppy);
				bgq_weylsite_double *weylsite = (bgq_weylsite_double *)((uint8_t*)&field->sec_collapsed_double[ic] - 32);

				#define PRECISION double
				#define BGQ_READWEYLLAYOUT_INC_
				#include "bgq_ReadWeyllayout.inc.c"
				#undef PRECISION
			}
		} else {
			//assert(field->hasFullspinorData);
			if (sloppy) {
				bgq_su3_spinor_load_float(spinor, &field->sec_fullspinor_float[ic]);
						bgq_spinorqpx_expect(spinor, t1, t2, x, y, z);
			} else {
				bgq_su3_spinor_load_double(spinor, &field->sec_fullspinor_double[ic]);
						bgq_spinorqpx_expect(spinor, t1, t2, x, y, z);
			}
		}

#if 0
		if (mul) {
			bgq_cmul(spinor_v0, qz, spinor_v0);
			bgq_cmul(spinor_v1, qz, spinor_v1);
			bgq_cmul(spinor_v2, qw, spinor_v2);
			bgq_cmul(spinor_v3, qw, spinor_v3);
		}
#endif
	}

	bgq_su3_spinor_mov(*target, spinor);
}


void bgq_spinorfield_prepareWrite(bgq_weylfield_controlblock *field, tristate isOdd, bgq_spinorfield_layout layout, bool preserveData);
void bgq_spinorfield_prepareReadWrite(bgq_weylfield_controlblock *field, tristate isOdd, bgq_spinorfield_layout layout);
bgq_spinorfield_layout bgq_spinorfield_prepareRead(bgq_weylfield_controlblock *field, tristate isOdd, bool acceptWeyl, bool acceptDouble, bool acceptFloat, bool acceptMul, bool acceptLegacy);

bgq_weylfield_controlblock *bgq_translate_spinorfield(const spinor *legacyField);
//void bgq_spinorfield_prepareLegacy(bgq_weylfield_controlblock *field, bool read);

void spinorfield_enable(const spinor *legacyField, int read, int write);
void spinorfield_propagateOddness(const spinor *targetLegacyField, const spinor *sourceLegacyField);
void spinorfield_propagateInvertedOddness(const spinor *targetLegacyField, const spinor *sourceLegacyField);


#ifdef BGQ_COORDCHECK
#define bgq_legacy_markcoords(isOdd, legacyField) bgq_legacy_markcoords_raw(isOdd, legacyField);
#else
#define bgq_legacy_markcoords(isOdd, legacyField)
#endif
void bgq_legacy_markcoords_raw(bool isOdd, spinor *legacyField);

void bgq_spinorfield_zero(bgq_weylfield_controlblock *field, tristate isOdd);


EXTERN_INLINE bool bgq_spinorfield_isOdd(bgq_weylfield_controlblock *field) {
	assert(field);
	if (field->has_fulllayout_double || field->has_fulllayout_float || field->has_weyllayout_double || field->has_weyllayout_float) {
		// field->isOdd is well-defined
	} else if (field->has_legacy) {
		assert(false);
	} else {
		assert(false);
	}
	return field->isOdd;
}

void bgq_spinorfield_annotateOddness(bgq_weylfield_controlblock *field, bool isOdd);
void spinorfield_setOddness(const spinor *field, int isOdd);
void bgq_spinorfield_dump(bgq_weylfield_controlblock *field, char *desc);
void spinorfield_dump(const spinor *field, char *desc);


#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_SPINORFIELD_H_ */
