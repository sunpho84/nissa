/*
 * bgq_gaugefield.h
 *
 *  Created on: Oct 26, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_GAUGEFIELD_H_
#define BGQ_GAUGEFIELD_H_

#include "bgq_utils.h"
#include "bgq_field.h"

#ifndef BGQ_GAUGEFIELD_C_
#define EXTERN_INLINE EXTERN_INLINE_DECLARATION
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE EXTERN_INLINE_DEFINITION
#define EXTERN_FIELD
#define EXTERN_INIT(val) = (val)
#endif

typedef struct {
	COMPLEX_PRECISION c[3][3];
} bgq_su3matrix;

typedef struct {
	COMPLEX_PRECISION c[3][3][PHYSICAL_LK];
} bgq_su3matrix_vec;

typedef struct {
	COMPLEX_PRECISION c[3][3][PHYSICAL_LK]; /* 3*3*2*sizeof(COMPLEX_PRECISION) = 288;144 bytes (4.5;2.25 L1 cache lines) */
	//COMPLEX_PRECISION padding[6];
} bgq_gaugesu3;
typedef struct {
	bgq_gaugesu3 su3[8];
} bgq_gaugesite;
typedef struct {
	bgq_gaugesu3 *(eodir[PHYSICAL_LP][PHYSICAL_LD]);
} bgq_gaugeeodir;
typedef bgq_gaugeeodir (*bgq_gaugefield);

EXTERN_FIELD bgq_gaugesite *g_bgq_gaugefield_fromCollapsed[PHYSICAL_LP];
EXTERN_FIELD bgq_gaugesite *g_bgq_gaugefield_fromHalfvolume[PHYSICAL_LP]; //deprecated //TODO: Remove, not needed
EXTERN_FIELD bgq_gaugesite *g_bgq_gaugefield_fromSurface[PHYSICAL_LP];//deprecated
EXTERN_FIELD bgq_gaugesite *g_bgq_gaugefield_fromBody[PHYSICAL_LP];//deprecated

#ifdef __cplusplus
extern "C" void bgq_gaugefield_init();
extern "C" void bgq_gaugefield_transferfrom(su3 *sourcefield);
#else
void bgq_gaugefield_init();
void bgq_gaugefield_transferfrom(su3 *sourcefield);
#endif

#define bgq_gauge_fromqpx(gaugeqpx,k) bgq_gauge_fromqpx_raw(bgq_su3_mvars(gaugeqpx),k)
EXTERN_INLINE bgq_su3matrix bgq_gauge_fromqpx_raw(bgq_su3_mparams(gaugeqpx), ucoord k) {
	assert(0 <= k && k < PHYSICAL_LK);
	bgq_su3matrix result = {{{0}}};
	result.c[0][0] = bgq_cmplxval(gaugeqpx_c00,k);
	result.c[0][1] = bgq_cmplxval(gaugeqpx_c01,k);
	result.c[0][2] = bgq_cmplxval(gaugeqpx_c02,k);
	result.c[1][0] = bgq_cmplxval(gaugeqpx_c10,k);
	result.c[1][1] = bgq_cmplxval(gaugeqpx_c11,k);
	result.c[1][2] = bgq_cmplxval(gaugeqpx_c12,k);
	result.c[2][0] = bgq_cmplxval(gaugeqpx_c20,k);
	result.c[2][1] = bgq_cmplxval(gaugeqpx_c21,k);
	result.c[2][2] = bgq_cmplxval(gaugeqpx_c22,k);
	return result;
}

void bgq_gauge_expect(bgq_su3matrix gauge,ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc);

#ifdef BGQ_COORDCHECK
#define bgq_gaugeqpx_expect(gaugeqpx,t1,t2,x,y,z,d,isSrc) bgq_gaugeqpx_expect_raw(bgq_su3_mvars(gaugeqpx),t1,t2,x,y,z,d,isSrc)
#else
#define bgq_gaugeqpx_expect(gaugeqpx,t1,t2,x,y,z,d,isSrc)
#endif
EXTERN_INLINE void bgq_gaugeqpx_expect_raw(bgq_su3_mparams(gaugeqpx),ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z, bgq_direction d,bool isSrc) {
	bgq_gauge_expect(bgq_gauge_fromqpx(gaugeqpx,0),t1,x,y,z,d,isSrc);
	bgq_gauge_expect(bgq_gauge_fromqpx(gaugeqpx,1),t2,x,y,z,d,isSrc);
}





#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_GAUGEFIELD_H_ */
