/*
 * bgq_qpx.c
 *
 *  Created on: Oct 26, 2012
 *      Author: meinersbur
 */

#define BGQ_QPX_C_
#include "bgq_qpx.h"

#include <stdbool.h>
#include <stdint.h>
#include <inttypes.h>

#define QPX_EXPECT(q, r0, i0, r1, i1) \
	do { \
		if (((int)bgq_elem0(q) != (int)(r0)) || ((int)bgq_elem1(q) != (int)(i0)) || ((int)bgq_elem2(q) != (int)(r1)) || ((int)bgq_elem3(q) != (int)(i1))) {  \
			failed = true; \
			fprintf(stderr, "QPX fail: (%g,%g,%g,%g) found, (%d,%d,%d,%d) expected at %s:%d\n",bgq_elem0(q),bgq_elem1(q),bgq_elem2(q),bgq_elem3(q),(int)(r0),(int)(i0),(int)(r1),(int)(i1), __FILE__,__LINE__); \
		} \
	} while (0)

#define VAL_EXPECT(val, exp) \
	do { \
		if ((val) != (exp)) {  \
			failed = true; \
			fprintf(stderr, "VAL fail: (%ld) found, (%ld) expected at %s:%d\n",(int64_t)(val),(int64_t)(exp), __FILE__,__LINE__); \
		} \
	} while (0)

void bgq_qpx_unittest(void) {
	printf("Running QPX unittest...\n");

	bool failed = false;
	complex_double data[2*4*3] __attribute__((aligned(32))) = {
			1 + 2*_Complex_I,
			3 + 4*_Complex_I,
			5 + 6*_Complex_I,
			7 + 8*_Complex_I,
			9 + 10*_Complex_I,
			11 + 12*_Complex_I,
			13 + 14*_Complex_I,
			15 + 16*_Complex_I,
			17 + 18*_Complex_I,
			19 + 20*_Complex_I,
			21 + 22*_Complex_I,
			23 + 24*_Complex_I,
			25 + 26*_Complex_I,
			27 + 28*_Complex_I,
			29 + 30*_Complex_I,
			31 + 32*_Complex_I,
			33 + 34*_Complex_I,
			35 + 36*_Complex_I,
			37 + 38*_Complex_I,
			39 + 40*_Complex_I,
			41 + 42*_Complex_I,
			43 + 44*_Complex_I,
			45 + 46*_Complex_I,
			47 + 48*_Complex_I
	};
	if (((uintptr_t)&data[0] & (32-1)) != 0) {
		failed = true;
		fprintf(stderr, "QPX fail: not 32-byte aligned at %s:%d\n", __FILE__, __LINE__);
	}

	if ((int)creal(data[0])!=(int)1) {
		failed = true;
		fprintf(stderr, "QPX fail: (%d) found, (%d) expected at %s:%d\n", (int)(creal(data[0])), 1, __FILE__, __LINE__);
	}




	bgq_vector4double_decl(v0);
	bgq_lda_double(v0, 0, &data[0]);
	QPX_EXPECT(v0,1,2,3,4);

	bgq_vector4double_decl(v1);
	bgq_lda_double(v1, 32, &data[0]);
	QPX_EXPECT(v1,5,6,7,8);

	bgq_vector4double_decl(v3);
	bgq_ld2a_double(v3, 0, &data[0]);
	QPX_EXPECT(v3,1,2,1,2);

	bgq_vector4double_decl(v4);
	bgq_ld2a_double(v4, 16, &data[0]);
	QPX_EXPECT(v4,3,4,3,4);

	bgq_vector4double_decl(lmerge);
	bgq_lmerge(lmerge, v0, v1);
	QPX_EXPECT(lmerge,1,2,5,6);

	bgq_vector4double_decl(rmerge);
	bgq_rmerge(rmerge, v0, v1);
	QPX_EXPECT(rmerge,3,4,7,8);

	bgq_vector4double_decl(vmerge2);
	bgq_merge2(vmerge2, v0, v1);
	QPX_EXPECT(vmerge2,3,4,5,6);

	bgq_vector4double_decl(vsplat);
	bgq_complxval_splat(vsplat, 1 + 2*_Complex_I);
	QPX_EXPECT(vsplat,1,2,1,2);

	{
		bgq_vector4double_decl(qvlfdxa);
		uint8_t * const addr = (void*)&data[0];
		bgq_qvlfdxa(qvlfdxa, addr, 32);
		QPX_EXPECT(qvlfdxa,5,6,7,8);
	}

	{
		bgq_vector4double_decl(qvlfduxa);
		uint8_t *addr = (void*)&data[0];
		bgq_qvlfduxa(qvlfduxa, addr, 32);
		QPX_EXPECT(qvlfduxa,5,6,7,8);
		if ((void*)addr != (void*)&data[2]) {
			failed = true;
			fprintf(stderr, "QPX fail: (%lu) found, (%lu) expected at %s:%d\n", (uintptr_t)addr, (uintptr_t)(&data[2]), __FILE__, __LINE__);
		}
	}

	bgq_su3_weyl_decl(weyl);
	bgq_su3_weyl_load_double(weyl, &data[0]);
	QPX_EXPECT(weyl_v0_c0,1,2,3,4);
	QPX_EXPECT(weyl_v0_c1,5,6,7,8);
	QPX_EXPECT(weyl_v0_c2,9,10,11,12);
	QPX_EXPECT(weyl_v1_c0,13,14,15,16);
	QPX_EXPECT(weyl_v1_c1,17,18,19,20);
	QPX_EXPECT(weyl_v1_c2,21,22,23,24);

	bgq_su3_spinor_decl(spinor);
	bgq_su3_spinor_load_double(spinor, &data[0]);
	QPX_EXPECT(spinor_v0_c0,1,2,3,4);
	QPX_EXPECT(spinor_v0_c1,5,6,7,8);
	QPX_EXPECT(spinor_v0_c2,9,10,11,12);
	QPX_EXPECT(spinor_v1_c0,13,14,15,16);
	QPX_EXPECT(spinor_v1_c1,17,18,19,20);
	QPX_EXPECT(spinor_v1_c2,21,22,23,24);
	QPX_EXPECT(spinor_v2_c0,25,26,27,28);
	QPX_EXPECT(spinor_v2_c1,29,30,31,32);
	QPX_EXPECT(spinor_v2_c2,33,34,35,36);
	QPX_EXPECT(spinor_v3_c0,37,38,39,40);
	QPX_EXPECT(spinor_v3_c1,41,42,43,44);
	QPX_EXPECT(spinor_v3_c2,45,46,47,48);

	double store[2*24] __attribute__((aligned(32))) = {0};
	if ((((uintptr_t)&store[0]) & (32-1)) != 0) {
		failed = true;
		fprintf(stderr, "QPX fail: not 32-byte aligned at %s:%d\n", __FILE__, __LINE__);
	}

	bgq_su3_weyl_store_double(&store, weyl);
	for (unsigned i = 0; i < 2*3*4;i+=1) {
		if ((int)store[i] != (int)(i+1)) {
			failed = true;
			fprintf(stderr, "QPX fail i=%d: (%d) found, (%d) expected at %s:%d\n", i, (int)(store[i]), (int)(i+1), __FILE__, __LINE__);
		}
	}

	bgq_su3_mdecl(matrix);
	bgq_su3_matrix_load_double(matrix, &data[0]);
	QPX_EXPECT(matrix_c00,1,2,3,4);
	QPX_EXPECT(matrix_c01,5,6,7,8);
	QPX_EXPECT(matrix_c02,9,10,11,12);
	QPX_EXPECT(matrix_c10,13,14,15,16);
	QPX_EXPECT(matrix_c11,17,18,19,20);
	QPX_EXPECT(matrix_c12,21,22,23,24);
	QPX_EXPECT(matrix_c20,25,26,27,28);
	QPX_EXPECT(matrix_c21,29,30,31,32);
	QPX_EXPECT(matrix_c22,33,34,35,36);

	{
		bgq_vector4double_decl(v_qvlfcduxa);
		bgq_zero(v_qvlfcduxa);
		uint8_t *addr = (void*)&data[0];
		bgq_qvlfcduxa(v_qvlfcduxa, addr, 32);
		QPX_EXPECT(v_qvlfcduxa,5,6,5,6);
		VAL_EXPECT(addr, (void*)&data[2]);
	}

	{
		uint8_t *addr = (void*)&store[0];
		bgq_qvstfcduxa(rmerge,addr,32);
		VAL_EXPECT(store[4], 3);
		VAL_EXPECT(store[5], 4);
		VAL_EXPECT(store[6], 7);
		VAL_EXPECT(store[7], 8);
		VAL_EXPECT(addr, (void*)&store[4]);
	}



// float precision ////////////////////////////////////////////////////////////

	complex_float data_float[2*4*3] __attribute__((aligned(16))) = {
			1 + 2*_Complex_I,
			3 + 4*_Complex_I,
			5 + 6*_Complex_I,
			7 + 8*_Complex_I,
			9 + 10*_Complex_I,
			11 + 12*_Complex_I,
			13 + 14*_Complex_I,
			15 + 16*_Complex_I,
			17 + 18*_Complex_I,
			19 + 20*_Complex_I,
			21 + 22*_Complex_I,
			23 + 24*_Complex_I,
			25 + 26*_Complex_I,
			27 + 28*_Complex_I,
			29 + 30*_Complex_I,
			31 + 32*_Complex_I,
			33 + 34*_Complex_I,
			35 + 36*_Complex_I,
			37 + 38*_Complex_I,
			39 + 40*_Complex_I,
			41 + 42*_Complex_I,
			43 + 44*_Complex_I,
			45 + 46*_Complex_I,
			47 + 48*_Complex_I
	};
	if (((uintptr_t)&data_float[0] & (16-1)) != 0) {
		failed = true;
		fprintf(stderr, "QPX fail: not 16-byte aligned at %s:%d\n", __FILE__, __LINE__);
	}


	bgq_vector4double_decl(v_qvlfsxa);
	bgq_qvlfsxa(v_qvlfsxa, &data_float[0], 0);
	QPX_EXPECT(v_qvlfsxa,1,2,3,4);

	{
		bgq_vector4double_decl(v_qvlfsuxa);
		uint8_t *addr = (uint8_t*)&data_float[0];
		bgq_qvlfsuxa(v_qvlfsuxa, addr, 16);
		QPX_EXPECT(v_qvlfsuxa,5,6,7,8);
		VAL_EXPECT(addr, (uint8_t*)&data_float[2]);
	}


	float store_float[2*24] __attribute__((aligned(16))) = {0};
	if ((((uintptr_t)&store_float[0]) & (16-1)) != 0) {
		failed = true;
		fprintf(stderr, "QPX fail: not 16-byte aligned at %s:%d\n", __FILE__, __LINE__);
	}

	bgq_qvstfsxa(v1, &store_float[0], 0);
	VAL_EXPECT(store_float[0], 5);
	VAL_EXPECT(store_float[1], 6);
	VAL_EXPECT(store_float[2], 7);
	VAL_EXPECT(store_float[3], 8);

	memset(&store_float, 0, sizeof(store_float));
	{
		uint8_t *addr = (uint8_t*)&store_float[0];
		bgq_qvstfsuxa(v0, addr, 16);
		VAL_EXPECT(store_float[4], 1);
		VAL_EXPECT(store_float[5], 2);
		VAL_EXPECT(store_float[6], 3);
		VAL_EXPECT(store_float[7], 4);
		VAL_EXPECT(addr,(uint8_t*)&store_float[4]);
	}


	bgq_su3_weyl_load_float(weyl, &data_float[0]);
	QPX_EXPECT(weyl_v0_c0,1,2,3,4);
	QPX_EXPECT(weyl_v0_c1,5,6,7,8);
	QPX_EXPECT(weyl_v0_c2,9,10,11,12);
	QPX_EXPECT(weyl_v1_c0,13,14,15,16);
	QPX_EXPECT(weyl_v1_c1,17,18,19,20);
	QPX_EXPECT(weyl_v1_c2,21,22,23,24);

	bgq_su3_spinor_load_float(spinor, &data_float[0]);
	QPX_EXPECT(spinor_v0_c0,1,2,3,4);
	QPX_EXPECT(spinor_v0_c1,5,6,7,8);
	QPX_EXPECT(spinor_v0_c2,9,10,11,12);
	QPX_EXPECT(spinor_v1_c0,13,14,15,16);
	QPX_EXPECT(spinor_v1_c1,17,18,19,20);
	QPX_EXPECT(spinor_v1_c2,21,22,23,24);
	QPX_EXPECT(spinor_v2_c0,25,26,27,28);
	QPX_EXPECT(spinor_v2_c1,29,30,31,32);
	QPX_EXPECT(spinor_v2_c2,33,34,35,36);
	QPX_EXPECT(spinor_v3_c0,37,38,39,40);
	QPX_EXPECT(spinor_v3_c1,41,42,43,44);
	QPX_EXPECT(spinor_v3_c2,45,46,47,48);

	memset(&store_float[0], 0, sizeof(store_float));
	bgq_su3_weyl_store_float(&store_float[0], weyl);
	for (unsigned i = 0; i < 2*3*4;i+=1) {
		if (store_float[i] != (i+1)) {
			failed = true;
			fprintf(stderr, "QPX fail i=%d: (%g) found, (%d) expected at %s:%d\n", i, store_float[i], i+1, __FILE__, __LINE__);
		}
	}

	memset(&store_float[0], 0, sizeof(store_float));
	bgq_su3_spinor_store_float(&store_float[0], spinor);
	for (unsigned i = 0; i < 4*3*4;i+=1) {
		if (store_float[i] != (i+1)) {
			failed = true;
			fprintf(stderr, "QPX fail i=%d: (%g) found, (%d) expected at %s:%d\n", i, store_float[i], i+1, __FILE__, __LINE__);
		}
	}


	{
		bgq_vector4double_decl(v_qvlfcsuxa);
		uint8_t *addr = (void*)&data_float[0];
		bgq_qvlfcsuxa(v_qvlfcsuxa, addr, 16);
		QPX_EXPECT(v_qvlfcsuxa,5,6,5,6);
		VAL_EXPECT(addr, (void*)&data_float[2]);
	}

	{
		uint8_t *addr = (void*)&store_float[0];
		bgq_qvstfcsuxa(rmerge,addr,16);
		VAL_EXPECT(store_float[4], 3);
		VAL_EXPECT(store_float[5], 4);
		VAL_EXPECT(store_float[6], 7);
		VAL_EXPECT(store_float[7], 8);
		VAL_EXPECT(addr, (void*)&store_float[4]);
	}

	/////////////////////////////////////////////////////////////////////////////


	if (!failed) {
		printf("QPX unittest passed\n");
	} else {
		//abort();
	}
}
