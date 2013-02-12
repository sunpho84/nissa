/*
 * bgq_qpx.h
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_QPX_H_
#define BGQ_QPX_H_

#include "bgq_utils.h"
//#include <mpi.h>
//#include "complex_c99.h"
#include <string.h>
#include <complex.h>
#include <stdint.h>

#if BGQ_QPX
#include <hwi/include/bqc/A2_inlines.h>
#else
#include <mpi.h>
#endif

#ifndef BGQ_QPX_C_
#define EXTERN_INLINE EXTERN_INLINE_DECLARATION
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE EXTERN_INLINE_DEFINITION
#define EXTERN_FIELD
#define EXTERN_INIT(val) = (val)
#endif







#ifndef BGQ_QPX
#define BGQ_QPX 0
#endif

#define BGQ_ALIGNMENT_QPX 32
#define BGQ_ALIGNMENT_L1 64
#define BGQ_ALIGNMENT_L1P 128
#define BGQ_ALIGNMENT_L2 128



typedef struct {
	double q[4];
} v4d;

typedef struct {
	float q[4];
} v4f;

typedef struct {
	_Complex double q[2];
} v2cd;

typedef struct {
	_Complex float q[2];
} v2cf;

#if defined(__IBMC__)||defined(__IBMCPP__)
	#define BGQ_VECTOR4DOUBLE_SUBSCRIPT(addr,idx) ( (*((vector4double*)(addr)))[(idx)] ) /* allows subscript like an array */
#else
	typedef v4d vector4double;
	#define BGQ_VECTOR4DOUBLE_SUBSCRIPT(addr,idx) ( ((v4d*)(addr))->q[(idx)] ) /* emulation using a struct */
#endif
#define BGQ_VECTOR4FLOAT_SUBSCRIPT(addr,idx) ( ((v4f*)(addr))->q[(idx)] ) /* emulation using a struct */




#if !BGQ_QPX
#define bgq_elem0(arg) \
	NAME2(arg,q0)

#define bgq_elem1(arg) \
	NAME2(arg,q1)

#define bgq_elem2(arg) \
	NAME2(arg,q2)

#define bgq_elem3(arg) \
	NAME2(arg,q3)

#define bgq_vars(name) \
	NAME2(name,q0),NAME2(name,q1),NAME2(name,q2),NAME2(name,q3)

#define bgq_types \
	double,double,double,double

#define bgq_params(name) \
	double NAME2(name,q0),double NAME2(name,q1),double NAME2(name,q2),double NAME2(name,q3)

#define bgq_complxval_splat(dst,cmplx) \
	bgq_elem0(dst) = creal(cmplx); \
	bgq_elem1(dst) = cimag(cmplx); \
	bgq_elem2(dst) = creal(cmplx); \
	bgq_elem3(dst) = cimag(cmplx)

#define bgq_vector4double_decl(name) \
	double NAME2(name,q0); \
	double NAME2(name,q1); \
	double NAME2(name,q2); \
	double NAME2(name,q3)

#define bgq_vector4double_decl_leftonly(name) \
	double NAME2(name,q0); \
	double NAME2(name,q1)

#define bgq_vector4double_decl_rightonly(name) \
	double NAME2(name,q2); \
	double NAME2(name,q3)

#define bgq_lda_double(dst,offset,addr)                                 \
	assert( (((uintptr_t)addr) + offset) % 32 == 0);                \
	NAME2(dst,q0) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 0); \
	NAME2(dst,q1) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 1); \
	NAME2(dst,q2) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 2); \
	NAME2(dst,q3) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 3)
//TODO: setting the 5 least significant bits of addr+offset to zero

#define bgq_lda_float(dst,offset,addr)                                 \
	assert( (((uintptr_t)addr) + offset) % 16 == 0);                \
	NAME2(dst,q0) = BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 0); \
	NAME2(dst,q1) = BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 1); \
	NAME2(dst,q2) = BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 2); \
	NAME2(dst,q3) = BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 3)

#define bgq_ld2a_double(dst,offset,addr)                                      \
	assert( (((uintptr_t)(addr)) + (offset)) % 16 == 0);                          \
	NAME2(dst,q0) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 0); \
	NAME2(dst,q1) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 1); \
	NAME2(dst,q2) = NAME2(dst,q0);                                             \
	NAME2(dst,q3) = NAME2(dst,q1)

#define bgq_ld2a_leftonly_double(dst,offset,addr)                                      \
	assert( (((uintptr_t)(addr)) + (offset)) % 16 == 0);                          \
	NAME2(dst,q0) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 0); \
	NAME2(dst,q1) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 1)

#define bgq_ld2a_float(dst,offset,addr) \
	assert( (((uintptr_t)(addr)) + (offset)) % 8 == 0);                         \
	NAME2(dst,q0) = BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 0); \
	NAME2(dst,q1) = BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 1); \
	NAME2(dst,q2) = NAME2(dst,q0);                                           \
	NAME2(dst,q3) = NAME2(dst,q1)

#define bgq_sta_double(src,offset,addr) \
	assert( (((uintptr_t)(addr)) + (offset)) % 32 == 0);                         \
	BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 0) = NAME2(src,q0); \
	BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 1) = NAME2(src,q1); \
	BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 2) = NAME2(src,q2); \
	BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 3) = NAME2(src,q3)

#define bgq_sta_float(src,offset,addr) \
	assert( (((uintptr_t)(addr)) + (offset)) % 16 == 0);                        \
	BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 0) = NAME2(src,q0); \
	BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 1) = NAME2(src,q1); \
	BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 2) = NAME2(src,q2); \
	BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 3) = NAME2(src,q3)

#define bgq_st2a_double(src,offset,addr) \
	assert( (((uintptr_t)(addr)) + (offset)) % 16 == 0); \
	BGQ_VECTOR4DOUBLE_SUBSCRIPT((uint8_t*)(addr) + (size_t)(offset), 0) = NAME2(src,q0); \
	BGQ_VECTOR4DOUBLE_SUBSCRIPT((uint8_t*)(addr) + (size_t)(offset), 1) = NAME2(src,q1)

#define bgq_st2a_float(src,offset,addr) \
	assert( (((uintptr_t)(addr)) + (offset)) % 8 == 0);                        \
	BGQ_VECTOR4FLOAT_SUBSCRIPT((uint8_t*)(addr) + (offset), 0) = NAME2(src,q0); \
	BGQ_VECTOR4FLOAT_SUBSCRIPT((uint8_t*)(addr) + (offset), 1) = NAME2(src,q1)

#define bgq_neg(dst,arg)            \
	NAME2(dst,q0) = - NAME2(arg,q0); \
	NAME2(dst,q1) = - NAME2(arg,q1); \
	NAME2(dst,q2) = - NAME2(arg,q2); \
	NAME2(dst,q3) = - NAME2(arg,q3)

#define bgq_add(dst,lhs,rhs)        \
	NAME2(dst,q0) = NAME2(lhs,q0) + NAME2(rhs,q0); \
	NAME2(dst,q1) = NAME2(lhs,q1) + NAME2(rhs,q1); \
	NAME2(dst,q2) = NAME2(lhs,q2) + NAME2(rhs,q2); \
	NAME2(dst,q3) = NAME2(lhs,q3) + NAME2(rhs,q3)

#define bgq_sub(dst,lhs,rhs)        \
	NAME2(dst,q0) = NAME2(lhs,q0) - NAME2(rhs,q0); \
	NAME2(dst,q1) = NAME2(lhs,q1) - NAME2(rhs,q1); \
	NAME2(dst,q2) = NAME2(lhs,q2) - NAME2(rhs,q2); \
	NAME2(dst,q3) = NAME2(lhs,q3) - NAME2(rhs,q3)

#define bgq_xxnpmadd(dst,a,b,c)                     \
	{                                               \
		bgq_vector4double_decl(MAKENAME5(xxnpmadd,dst,a,b,v)); \
		MAKENAME6(xxnpmadd,dst,a,b,v,q0) = - (CONCAT(a,_q1) * CONCAT(b,_q1) - CONCAT(c,_q0)); \
		MAKENAME6(xxnpmadd,dst,a,b,v,q1) =    CONCAT(a,_q0) * CONCAT(b,_q1) + CONCAT(c,_q1) ; \
		MAKENAME6(xxnpmadd,dst,a,b,v,q2) = - (CONCAT(a,_q3) * CONCAT(b,_q3) - CONCAT(c,_q2)); \
		MAKENAME6(xxnpmadd,dst,a,b,v,q3) =    CONCAT(a,_q2) * CONCAT(b,_q3) + CONCAT(c,_q3) ; \
		bgq_mov(dst,MAKENAME5(xxnpmadd,dst,a,b,v));                                           \
	}

#define bgq_iadd(dst,lhs,rhs)         \
	{                                 \
		bgq_vector4double_decl(tmp);  \
		tmp_q0 = NAME2(lhs,q0) - NAME2(rhs,q1); \
		tmp_q1 = NAME2(lhs,q1) + NAME2(rhs,q0); \
		tmp_q2 = NAME2(lhs,q2) - NAME2(rhs,q3); \
		tmp_q3 = NAME2(lhs,q3) + NAME2(rhs,q2); \
		bgq_mov(dst,tmp);             \
	}

#define bgq_isub(dst,lhs,rhs)         \
	{                                 \
		bgq_vector4double_decl(tmp);  \
		tmp_q0 = NAME2(lhs,q0) + NAME2(rhs,q1); \
		tmp_q1 = NAME2(lhs,q1) - NAME2(rhs,q0); \
		tmp_q2 = NAME2(lhs,q2) + NAME2(rhs,q3); \
		tmp_q3 = NAME2(lhs,q3) - NAME2(rhs,q2); \
		bgq_mov(dst,tmp);             \
	}

#define bgq_mul(dst,lhs,rhs)                                          \
	{                                                                  \
		bgq_vector4double_decl(MAKENAME4(mul,dst,lhs,rhs));            \
		MAKENAME5(mul,dst,lhs,rhs,q0) = NAME2(lhs,q0) * NAME2(rhs,q0); \
		MAKENAME5(mul,dst,lhs,rhs,q1) = NAME2(lhs,q1) * NAME2(rhs,q1); \
		MAKENAME5(mul,dst,lhs,rhs,q2) = NAME2(lhs,q2) * NAME2(rhs,q2); \
		MAKENAME5(mul,dst,lhs,rhs,q3) = NAME2(lhs,q3) * NAME2(rhs,q3); \
		bgq_mov(dst,MAKENAME4(mul,dst,lhs,rhs));                       \
	}

#define bgq_merge2(dst,a23_to01,b01_to23) \
	{                                       \
		bgq_vector4double_decl(MAKENAME4(merge2,dst,a23_to01,b01_to23));      \
		MAKENAME5(merge2,dst,a23_to01,b01_to23,q0) = NAME2(a23_to01,q2);      \
		MAKENAME5(merge2,dst,a23_to01,b01_to23,q1) = NAME2(a23_to01,q3);      \
		MAKENAME5(merge2,dst,a23_to01,b01_to23,q2) = NAME2(b01_to23,q0);      \
		MAKENAME5(merge2,dst,a23_to01,b01_to23,q3) = NAME2(b01_to23,q1);      \
		bgq_mov(dst,MAKENAME4(merge2,dst,a23_to01,b01_to23));                   \
	}

// returns the left complex part of both arguments, i.e.
// (a_q0 a_q1 b_q0 b_q1)
#define bgq_lmerge(dst, a01_to01, b01_to23) \
	do {                                     \
		bgq_vector4double_decl(tmp);         \
		NAME2(tmp,q0) = NAME2(a01_to01,q0);  \
		NAME2(tmp,q1) = NAME2(a01_to01,q1);  \
		NAME2(tmp,q2) = NAME2(b01_to23,q0);  \
		NAME2(tmp,q3) = NAME2(b01_to23,q1);  \
		bgq_mov(dst,tmp);                    \
	} while (0)

// returns the right complex part of both arguments, i.e.
// (a_q2 a_q3 b_q2 b_q3)
#define bgq_rmerge(dst, a23_to01, b23_to23) \
	do {                                     \
		bgq_vector4double_decl(tmp);         \
		NAME2(tmp,q0) = NAME2(a23_to01,q2);  \
		NAME2(tmp,q1) = NAME2(a23_to01,q3);  \
		NAME2(tmp,q2) = NAME2(b23_to23,q2);  \
		NAME2(tmp,q3) = NAME2(b23_to23,q3);  \
		bgq_mov(dst,tmp);                    \
	} while (0)

#define bgq_xmul(dst,lhs,rhs)                                          \
	do {                                                                \
		bgq_vector4double_decl(MAKENAME4(xmul,dst,lhs,rhs));            \
		MAKENAME5(xmul,dst,lhs,rhs,q0) = NAME2(lhs,q0) * NAME2(rhs,q0); \
		MAKENAME5(xmul,dst,lhs,rhs,q1) = NAME2(lhs,q0) * NAME2(rhs,q1); \
		MAKENAME5(xmul,dst,lhs,rhs,q2) = NAME2(lhs,q2) * NAME2(rhs,q2); \
		MAKENAME5(xmul,dst,lhs,rhs,q3) = NAME2(lhs,q2) * NAME2(rhs,q3); \
		bgq_mov(dst,MAKENAME4(xmul,dst,lhs,rhs));                       \
	} while (0)

#define bgq_xmadd(dst,a,b,c)                                                    \
	do {                                                                         \
		bgq_vector4double_decl(MAKENAME5(xmadd,dst,a,b,c));                      \
		MAKENAME6(xmadd,dst,a,b,c,q0) = NAME2(a,q0) * NAME2(b,q0) + NAME2(c,q0); \
		MAKENAME6(xmadd,dst,a,b,c,q1) = NAME2(a,q0) * NAME2(b,q1) + NAME2(c,q1); \
		MAKENAME6(xmadd,dst,a,b,c,q2) = NAME2(a,q2) * NAME2(b,q2) + NAME2(c,q2); \
		MAKENAME6(xmadd,dst,a,b,c,q3) = NAME2(a,q2) * NAME2(b,q3) + NAME2(c,q3); \
		bgq_mov(dst,MAKENAME5(xmadd,dst,a,b,c));                                 \
	} while (0)

#define bgq_madd(dst,a,b,c)                                                    \
	{                                                                           \
		bgq_vector4double_decl(MAKENAME5(madd,dst,a,b,c));                      \
		MAKENAME6(madd,dst,a,b,c,q0) = NAME2(a,q0) * NAME2(b,q0) + NAME2(c,q0); \
		MAKENAME6(madd,dst,a,b,c,q1) = NAME2(a,q1) * NAME2(b,q1) + NAME2(c,q1); \
		MAKENAME6(madd,dst,a,b,c,q2) = NAME2(a,q2) * NAME2(b,q2) + NAME2(c,q2); \
		MAKENAME6(madd,dst,a,b,c,q3) = NAME2(a,q3) * NAME2(b,q3) + NAME2(c,q3); \
		bgq_mov(dst,MAKENAME5(madd,dst,a,b,c));                                 \
	}

#define bgq_msub(dst,a,b,c)                                                    \
	{                                                                           \
		bgq_vector4double_decl(MAKENAME5(madd,dst,a,b,c));                      \
		MAKENAME6(madd,dst,a,b,c,q0) = NAME2(a,q0) * NAME2(b,q0) - NAME2(c,q0); \
		MAKENAME6(madd,dst,a,b,c,q1) = NAME2(a,q1) * NAME2(b,q1) - NAME2(c,q1); \
		MAKENAME6(madd,dst,a,b,c,q2) = NAME2(a,q2) * NAME2(b,q2) - NAME2(c,q2); \
		MAKENAME6(madd,dst,a,b,c,q3) = NAME2(a,q3) * NAME2(b,q3) - NAME2(c,q3); \
		bgq_mov(dst,MAKENAME5(madd,dst,a,b,c));                                 \
	}

#define bgq_mov(dst,src) \
	NAME2(dst,q0) = NAME2(src,q0); \
	NAME2(dst,q1) = NAME2(src,q1); \
	NAME2(dst,q2) = NAME2(src,q2); \
	NAME2(dst,q3) = NAME2(src,q3)

#define bgq_xxcpnmadd(dst,a,b,c)                      \
	{                                                 \
		bgq_vector4double_decl(MAKENAME5(xmul,dst,a,b,c));   \
		MAKENAME6(xmul,dst,a,b,c,q0) =    CONCAT(a,_q1) * CONCAT(b,_q1) + CONCAT(c,_q0) ; \
		MAKENAME6(xmul,dst,a,b,c,q1) = - (CONCAT(a,_q0) * CONCAT(b,_q1) - CONCAT(c,_q1)); \
		MAKENAME6(xmul,dst,a,b,c,q2) =    CONCAT(a,_q3) * CONCAT(b,_q3) + CONCAT(c,_q2) ; \
		MAKENAME6(xmul,dst,a,b,c,q3) = - (CONCAT(a,_q2) * CONCAT(b,_q3) - CONCAT(c,_q3)); \
		bgq_mov(dst,MAKENAME5(xmul,dst,a,b,c));                                           \
	}

#define bgq_cconst(dst,re,im) \
	NAME2(dst,q0) = re;            \
	NAME2(dst,q1) = im;            \
	NAME2(dst,q2) = re;            \
	NAME2(dst,q3) = im

#define bgq_zero(dst) \
	NAME2(dst,q0) = 0;     \
	NAME2(dst,q1) = 0;     \
	NAME2(dst,q2) = 0;     \
	NAME2(dst,q3) = 0


#define bgq_qvlfdxa(dst,addr,offset) \
	bgq_lda_double(dst,offset,addr)
#define bgq_qvlfduxa(dst,addr,offset) \
	do { \
		(addr) = (void*)((uintptr_t)(addr) + (offset)); \
		bgq_lda_double(dst,0,addr); \
	} while (0)

#define bgq_qvlfsxa(dst,addr,offset) \
	bgq_lda_float(dst,offset,addr)
#define bgq_qvlfsuxa(dst,addr,offset) \
	do {                               \
		(addr) = (void*)((uintptr_t)(addr) + (offset)); \
		bgq_lda_float(dst,0,addr);    \
	} while (0)






#define bgq_qvstfdxa(data,addr,offset) \
	bgq_sta_double(data,offset,addr)
#define bgq_qvstfduxa(data,addr,offset) \
	(addr) = (void*)((uintptr_t)(addr) + (offset)); \
	bgq_sta_double(data,0,addr)

#define bgq_qvstfsxa(dst,addr,offset) \
	bgq_sta_float(dst,offset,addr)
#define bgq_qvstfsuxa(dst,addr,offset) \
	do {                               \
		(addr) = (void*)((uintptr_t)(addr) + (offset)); \
		bgq_sta_float(dst,0,addr);    \
	} while (0)



#define bgq_qvlfcsuxa(dst,addr,offset) \
	do { \
		(addr) = (void*)((char*)(addr) + (offset)); \
		bgq_ld2a_float(dst,0,addr); \
	} while (0)

#define bgq_qvlfcduxa(dst,addr,offset) \
	do { \
		(addr) = (void*)((char*)(addr) + (offset)); \
		bgq_ld2a_double(dst,0,addr); \
	} while (0)

#define bgq_qvstfcduxa(data,addr,offset) \
	(addr) = (void*)((char*)(addr) + (offset)); \
	bgq_st2a_double(data,0,addr)

#define bgq_qvstfcsuxa(data,addr,offset) \
	(addr) = (void*)((char*)(addr) + (offset)); \
	bgq_st2a_float(data,0,addr)

#define bgq_dcbz(addr,offset) \
	do { \
		assert(((((uintptr_t)(addr)+(offset)))&(BGQ_ALIGNMENT_L1-1)) == 0); \
		memset((void*)(((uintptr_t)(addr) + (offset) + BGQ_ALIGNMENT_L1-1)&~(BGQ_ALIGNMENT_L1-1)), 0x00, 64); \
	} while (0)


#define bgq_dcbt(addr,offset)
#define bgq_dcbt_0(addr)
#define bgq_dcbi(addr,offset)
#define bgq_dcbi_0(addr)

#else
#define bgq_elem0(arg) \
	((arg)[0])

#define bgq_elem1(arg) \
	((arg)[1])

#define bgq_elem2(arg) \
	((arg)[2])

#define bgq_elem3(arg) \
	((arg)[3])

#define bgq_vars(name) \
	name

#define bgq_params(name) \
	vector4double (name)

#define bgq_types \
		vector4double

#define bgq_vector4double_decl(name) \
	vector4double name

#define bgq_complxval_splat(dst,cmplx) \
	(dst) = (vector4double){creal(cmplx),cimag(cmplx),creal(cmplx),cimag(cmplx)}

#define bgq_lda_double(dst,offset,addr) \
	assert((((uintptr_t)(addr)+(uintptr_t)(offset)) & 31) == 0 && "Must be aligned to 32 bytes"); \
	(dst) = vec_lda(offset,(double*)(addr))

#define bgq_lda_float(dst,offset,addr) \
	(dst) = vec_lda(offset,(float*)(addr))

#define bgq_ld2a_double(dst,offset,addr) \
	assert((((uintptr_t)(addr)+(uintptr_t)(offset)) & 15) == 0 && "Must be aligned to 16 bytes"); \
	(dst) = vec_ld2a(offset, (double*)(addr))

#define bgq_ld2a_float(dst,offset,addr) \
	(dst) = vec_ld2a(offset, (float*)(addr))

#define bgq_sta_double(src,offset,addr) \
	vec_sta(src, offset, (double*)(addr))

#define bgq_sta_float(src,offset,addr) \
	vec_sta(src, offset, (float*)(addr))

#define bgq_st2a_double(src,offset,addr) \
	vec_st2a(src, offset, (double*)(addr))

#define bgq_st2a_float(src,offset,addr) \
	vec_st2a(src, offset, (float*)(addr))

#define bgq_neg(dst,arg) \
	(dst) = vec_neg(arg);

#define bgq_add(dst,lhs,rhs) \
	(dst) = vec_add(lhs, rhs)

#define bgq_sub(dst,lhs,rhs) \
	(dst) = vec_sub(lhs, rhs)

#define bgq_mul(dst,lhs,rhs) \
	(dst) = vec_mul(lhs,rhs)

#define bgq_xxnpmadd(dst,a,b,c) \
	(dst) = vec_xxnpmadd(a,b,c)

#define bgq_iadd(dst,lhs,rhs) \
	bgq_xxnpmadd(dst,rhs,(vector4double)(1),lhs)

#define bgq_isub(dst,lhs,rhs) \
	bgq_xxcpnmadd(dst,rhs,(vector4double)(1),lhs)

// returns the middle part of both arguments concatenated, i.e.
// a_q0 a_q1 (a_q2 a_q3 b_q0 b_q1) b_q2 b_q3
#define bgq_merge2(dst, a23_to01, b01_to23) \
	asm ("qvaligni %[qrt],%[qra],%[qrb],2" : [qrt] "=v" (dst) : [qra] "v" (a23_to01), [qrb] "v" (b01_to23))
//	(dst) = vec_perm(a23_to01, b01_to23, vec_gpci(02345))

// returns the left complex part of both arguments, i.e.
// (a_q0 a_q1 b_q0 b_q1)
#define bgq_lmerge(dst, a01_to01, b01_to23) \
	(dst) = vec_perm(a01_to01, b01_to23, vec_gpci(00145))

// returns the right complex part of both arguments, i.e.
// (a_q2 a_q3 b_q2 b_q3)
#define bgq_rmerge(dst, a23_to01, b23_to23) \
	(dst) = vec_perm(a23_to01, b23_to23, vec_gpci(02367))

#define bgq_xmul(dst,lhs,rhs) \
	(dst) = vec_xmul(lhs,rhs)

#define bgq_xmadd(dst,a,b,c) \
	(dst) = vec_xmadd(a,b,c)

#define bgq_madd(dst,a,b,c) \
	(dst) = vec_madd(a,b,c)

#define bgq_msub(dst,a,b,c) \
	(dst) = vec_msub(a,b,c)

#define bgq_mov(dst,src) \
	(dst) = (src)

#define bgq_xxcpnmadd(dst,a,b,c) \
	(dst) = vec_xxcpnmadd(a,b,c)

#define bgq_cconst(dst,re,im) \
	(dst) = (vector4double){re,im,re,im}

#define bgq_zero(dst) \
	(dst) = (vector4double)(0)
//	dst = (vector4double){0,0,0,0}
//	dst = vec_logical(dst, dst, 0);

#define bgq_qvlfdxa(dst,addr,offset) \
	assert((((uintptr_t)(addr)+(uintptr_t)(offset)) & 31) == 0 && "Must be aligned to 32 bytes"); \
	asm ("qvlfdxa %[v4d],%[ptr],%[off]  \n" : [v4d] "=v" (dst) : [ptr]  "b" (addr),  [off] "r" (offset) ) /* no memory clobber, so pay attention! */
#define bgq_qvlfduxa(dst,addr,offset) \
	assert((((uintptr_t)(addr)+(uintptr_t)(offset)) & 31) == 0 && "Must be aligned to 32 bytes"); \
	asm ("qvlfduxa %[v4d],%[ptr],%[off]  \n" : [v4d] "=v" (dst), [ptr] "+b" (addr) : [off] "r" (offset) ) /* no memory clobber, so pay attention! */

#define bgq_qvlfsxa(dst,addr,offset) \
	assert((((uintptr_t)(addr)+(uintptr_t)(offset)) & 15) == 0 && "Must be aligned to 16 bytes"); \
	asm ("qvlfsxa %[v4d],%[ptr],%[off]  \n" : [v4d] "=v" (dst) : [ptr]  "b" (addr),  [off] "r" (offset) ) /* no memory clobber, so pay attention! */
#define bgq_qvlfsuxa(dst,addr,offset) \
	assert((((uintptr_t)(addr)+(uintptr_t)(offset)) & 15) == 0 && "Must be aligned to 16 bytes"); \
	asm ("qvlfsuxa %[v4d],%[ptr],%[off]  \n" : [v4d] "=v" (dst), [ptr] "+b" (addr) : [off] "r" (offset) ) /* no memory clobber, so pay attention! */

#define bgq_qvstfdxa(data,addr,offset) \
	asm ("qvstfdxa %[v4d],%[ptr],%[off]  \n" : : [ptr]  "b" (addr), [v4d] "v" (data), [off] "r" (offset) ) /* no memory clobber, so pay attention! (i.e. do not read from the memory location written here) */
#define bgq_qvstfduxa(data,addr,offset) \
	asm ("qvstfduxa %[v4d],%[ptr],%[off]  \n" : [ptr] "+b" (addr) : [v4d] "v" (data), [off] "r" (offset) ) /* no memory clobber, so pay attention! (i.e. do not read from the memory location written here) */

#define bgq_qvstfsxa(data,addr,offset) \
	asm ("qvstfsxa %[v4d],%[ptr],%[off]  \n" : : [ptr]  "b" (addr), [v4d] "v" (data), [off] "r" (offset) ) /* no memory clobber, so pay attention! (i.e. do not read from the memory location written here) */
#define bgq_qvstfsuxa(data,addr,offset) \
	asm ("qvstfsuxa %[v4d],%[ptr],%[off]  \n" : [ptr] "+b" (addr) : [v4d] "v" (data), [off] "r" (offset) ) /* no memory clobber, so pay attention! (i.e. do not read from the memory location written here) */


#define bgq_qvlfcsuxa(dst,addr,offset) \
	asm ("qvlfcsuxa %[v4d],%[ptr],%[off]  \n" : [v4d] "=v" (dst), [ptr] "+b" (addr) : [off] "r" (offset) ) /* no memory clobber, so pay attention! (i.e. do not read from the memory location written here) */

#define bgq_qvlfcduxa(dst,addr,offset) \
	asm ("qvlfcduxa %[v4d],%[ptr],%[off]  \n" : [v4d] "=v" (dst), [ptr] "+b" (addr) : [off] "r" (offset) ) /* no memory clobber, so pay attention! (i.e. do not read from the memory location written here) */


#define bgq_qvstfcduxa(data,addr,offset) \
	asm ("qvstfcduxa %[v4d],%[ptr],%[off]  \n" : [ptr] "+b" (addr) : [v4d] "v" (data), [off] "r" (offset) )

#define bgq_qvstfcsuxa(data,addr,offset) \
	asm ("qvstfcsuxa %[v4d],%[ptr],%[off]  \n" : [ptr] "+b" (addr) : [v4d] "v" (data), [off] "r" (offset) )


#define bgq_dcbt(addr,offset)  \
	asm ("dcbt %[offset],%[ptr]  \n" : : [offset] "b" (offset), [ptr] "r" (addr)) /* no memory clobber, hope that the compiler doesn't move this around too much */

#define bgq_dcbt_0(addr) \
	asm ("dcbt 0,%[ptr]  \n" : : [ptr] "r" (addr)) /* no memory clobber, hope that the compiler doesn't move this around too much */

#define bgq_dcbi(addr,offset)

// Causes SIGILL
#define bgq_dcbi(addr,offset)
	//asm ("dcbi %offset,%[ptr]  \n" : : [offset] "b" (offset), [ptr] "r" (addr)) /* no memory clobber, hope that the compiler doesn't move this around too much */

// Causes SIGILL
#define bgq_dcbi_0(addr)
	//asm ("dcbi 0,%[ptr]  \n" : : [ptr] "r" (addr)) /* no memory clobber, hope that the compiler doesn't move this around too much */

#define bgq_dcbz(addr,offset)  \
	asm ("dcbz %[offset],%[ptr]  \n" : : [offset] "b" (offset), [ptr] "r" (addr))

#define bgq_dcbz_0(addr)  \
	asm ("dcbz 0,%[ptr]  \n" : : [ptr] "r" (addr))

//#define bgq_antioptaway(qreg)
//	asm ("qvstfdxa %[v4d],0,%[ptr],  \n" : : [ptr] "r" (&dummi), [v4d] "v" (qreg) )
	//asm ("\n" : : [ptr] "v" (qreg))

//EXTERN_FIELD vector4double dummi;
//#define bgq_valgen(qreg)
//	asm ("qvgpci %[dst],0\n" : [dst] "=v" (qreg) : : )
	//asm ("\n" : : [ptr] "v" (qreg))

#endif


#define bgq_lda NAME2(bgq_lda,PRECISION)
#define bgq_sta NAME2(bgq_sta,PRECISION)
#define bgq_ld2a NAME2(bgq_ld2a,PRECISION)
#define bgq_st2a NAME2(bgq_st2a,PRECISION)

#define bgq_qvlfxa NAME2(bgq_qvlfxa,PRECISION)
#define bgq_qvlfxa_double bgq_qvlfdxa
#define bgq_qvlfxa_float  bgq_qvlfsxa

#define bgq_qvlfuxa NAME2(bgq_qvlfuxa,PRECISION)
#define bgq_qvlfuxa_double bgq_qvlfduxa
#define bgq_qvlfuxa_float bgq_qvlfsuxa

#define bgq_qvlfcuxa NAME2(bgq_qvlfcuxa,PRECISION)
#define bgq_qvlfcuxa_double bgq_qvlfcduxa
#define bgq_qvlfcuxa_float bgq_qvlfcsuxa

#define bgq_qvstfcuxa NAME2(bgq_qvstfcuxa,PRECISION)
#define bgq_qvstfcuxa_double bgq_qvstfcduxa
#define bgq_qvstfcuxa_float bgq_qvstfcsuxa

// vec_xmul(a, b)
// re =    a.re * b.re
// im =    a.re * b.im

// vec_xmul(b, a)
// re =    a.re * b.re
// im =    a.im * b.re

// vec_xmadd(a, b, c)
// re =    a.re * b.re + c.re
// im =    a.re * b.im + c.im

// vec_xmadd(b, a, c)
// re =    a.re * b.re + c.re
// im =    a.im * b.re + c.im

// vec_xxnpmadd(a, b, c)
// re =   -a.im * b.im + c.re
// im =    a.re * b.im + c.im

// vec_xxnpmadd(b, a, c)
// re =   -a.im * b.im + c.re
// im =    a.im * b.re + c.im

// vec_xxcpnmadd(a, b, c)
// re =    a.im * b.im + c.re
// im =   -a.re * b.im + c.im

// vec_xxcpnmadd(b, a, c)
#define bgq_rxxcpnmadd(dst, a, b, c) bgq_xxcpnmadd(dst, b, a, c)
// re =    a.im * b.im + c.re
// im =   -a.im * b.re + c.im

// vec_xxmadd(a, b, c)
// re =    a.im * b.im + c.re
// im =    a.re * b.im + c.im

// vec_xxmadd(b, a, c)
// re =    a.im * b.im + c.re
// im =    a.im * b.re + c.im



#define bgq_cmplxval1(name) \
	(bgq_elem0(name) + (bgq_elem1(name) * _Complex_I))

#define bgq_cmplxval2(name) \
	(bgq_elem2(name) + (bgq_elem3(name) * _Complex_I))

#define bgq_cmplxval(name,k) \
	((k) ? bgq_cmplxval2(name) : bgq_cmplxval1(name))





#define cvec_mul(a,b) vec_xxnpmadd(b,a,vec_xmul(a,b))
// vec_xxnpmadd(b,a,vec_xmul(a,b))
// vec_xxnpmadd(a,b,vec_xmul(b,a))
// re = a.re * b.re - a.im * b.im
// im = a.re * b.im + a.im * b.re
//    = a * b

#define bgq_cmul(dst,lhs,rhs)                                                                       \
	{                                                                                               \
		bgq_vector4double_decl(MAKENAME4(cmul,dst,lhs,rhs));                                         \
		bgq_xmul              (MAKENAME4(cmul,dst,lhs,rhs), lhs, rhs);                              \
		bgq_xxnpmadd          (dst                        , rhs, lhs, MAKENAME4(cmul,dst,lhs,rhs)); \
	}


// Conjugated lhs
// conj(lhs) * rhs
#define bgq_cjgmul(dst,lhs,rhs)                                                                       \
	do {                                                                                               \
		bgq_vector4double_decl(MAKENAME4(cjgmul,dst,lhs,rhs));                                         \
		bgq_xmul              (MAKENAME4(cjgmul,dst,lhs,rhs), lhs, rhs);                              \
		bgq_xxcpnmadd         (dst                          , rhs, lhs, MAKENAME4(cjgmul,dst,lhs,rhs)); \
	} while (0)


#define cvec_madd(a,b,c) vec_xxnpmadd(b,a,vec_xmadd(a,b,c))
// vec_xxnpmadd(b,a,vec_xmadd(a,b,c))
// vec_xxnpmadd(a,b,vec_xmadd(b,a,c))
// re = a.re * b.re - a.im * b.im + c.re
// im = a.re * b.im + a.im * b.re + c.im
//    = a * b + c

#define bgq_cmadd(dst,a,b,c)                                                                                  \
	{                                                                                                         \
		bgq_vector4double_decl(MAKENAME5(cmadd,dst,a,b,c));                                           \
		bgq_xmadd             (MAKENAME5(cmadd,dst,a,b,c), a, b, c);                                  \
		bgq_xxnpmadd          (dst                       , b, a, MAKENAME5(cmadd,dst,a,b,c)); \
	}

// conjugated second argument
#define cvec_pcmul(a,b) vec_xxcpnmadd(a,b,vec_xmul(b,a))
// vec_xxcpnmadd(a,b,vec_xmul(b,a))
// re = a.re * b.re + a.im * b.im
// im = a.im * b.re - a.re * b.im
//    = a * conjugate(b)

#define bgq_ccmul(dst,lhs,rhs) \
	{                                  \
		bgq_vector4double_decl(MAKENAME4(ccmul,dst,lhs,rhs));   \
		bgq_xmul(MAKENAME4(ccmul,dst,lhs,rhs),rhs,lhs);         \
		bgq_xxcpnmadd(dst,lhs,rhs,MAKENAME4(ccmul,dst,lhs,rhs)); \
	}

// conjugated second argument
#define cvec_pcpmadd(a,b,c) vec_xxcpnmadd(a,b,vec_xmadd(b,a,c))
// vec_xxcpnmadd(a,b,vec_xmadd(b,a,c))
// vec_xmadd(b,a,vec_xxcpnmadd(a,b,c))
// re = a.re * b.re + a.im * b.im + c.re
// im = a.im * b.re - a.re * b.im + c.im
//    = a * conjugate(b) + c

#define bgq_ccmadd(dst,a,b,c)        \
	{                                \
		bgq_vector4double_decl(MAKENAME5(ccmadd,dst,a,b,c)); \
		bgq_xmadd(MAKENAME5(ccmadd,dst,a,b,c),b,a,c);        \
		bgq_xxcpnmadd(dst,a,b,MAKENAME5(ccmadd,dst,a,b,c));  \
	}

#define cvec_piadd(a,b) vec_xxnpmadd((vector4double)(1),b,a)
// vec_xxnpmadd((vector4double)(1),b,a)
// re = a.re - b.im
// im = a.im + b.re
//    = a + bi

#define cvec_pisub(a,b) vec_xxcpnmadd(b,(vector4double)(1),a)
// vec_xxcpnmadd(b,(vector4double)(1),a)
// re = a.re + b.im
// im = a.im - b.re
//    = a - bi

// obsolete; use bgq_merge2
#define cvec_merge(a23_to01, b01_to23) vec_perm(a23_to01, b01_to23, vec_gpci(02345))
// [0] = a[2] = a.re
// [1] = a[3] = a.im
// [2] = b[0] = b.re
// [3] = b[1] = b.im

#define bgq_su3_vdecl(name)            \
	bgq_vector4double_decl(NAME2(name,c0)); \
	bgq_vector4double_decl(NAME2(name,c1)); \
	bgq_vector4double_decl(NAME2(name,c2))

#define bgq_su3_vdecl_leftonly(name)            \
	bgq_vector4double_decl_leftonly(NAME2(name,c0)); \
	bgq_vector4double_decl_leftonly(NAME2(name,c1)); \
	bgq_vector4double_decl_leftonly(NAME2(name,c2))

#define bgq_su3_vvars(name) \
	bgq_vars(NAME2(name,c0)), \
	bgq_vars(NAME2(name,c1)), \
	bgq_vars(NAME2(name,c2))

#define bgq_su3_vparams(name) \
	bgq_params(NAME2(name,c0)), \
	bgq_params(NAME2(name,c1)), \
	bgq_params(NAME2(name,c2))

#define bgq_su3_vrefparams(name) \
	bgq_refparams(NAME2(name,c0)), \
	bgq_refparams(NAME2(name,c1)), \
	bgq_refparams(NAME2(name,c2))

#define bgq_su3_mdecl(name)             \
	bgq_vector4double_decl(NAME2(name,c00)); \
	bgq_vector4double_decl(NAME2(name,c01)); \
	bgq_vector4double_decl(NAME2(name,c02)); \
	bgq_vector4double_decl(NAME2(name,c10)); \
	bgq_vector4double_decl(NAME2(name,c11)); \
	bgq_vector4double_decl(NAME2(name,c12)); \
	bgq_vector4double_decl(NAME2(name,c20)); \
	bgq_vector4double_decl(NAME2(name,c21)); \
	bgq_vector4double_decl(NAME2(name,c22))

#define bgq_su3_mvars(name) \
	bgq_vars(NAME2(name,c00)), \
	bgq_vars(NAME2(name,c01)), \
	bgq_vars(NAME2(name,c02)), \
	bgq_vars(NAME2(name,c10)), \
	bgq_vars(NAME2(name,c11)), \
	bgq_vars(NAME2(name,c12)), \
	bgq_vars(NAME2(name,c20)), \
	bgq_vars(NAME2(name,c21)), \
	bgq_vars(NAME2(name,c22))

#define bgq_su3_mparams(name) \
	bgq_params(NAME2(name,c00)), \
	bgq_params(NAME2(name,c01)), \
	bgq_params(NAME2(name,c02)), \
	bgq_params(NAME2(name,c10)), \
	bgq_params(NAME2(name,c11)), \
	bgq_params(NAME2(name,c12)), \
	bgq_params(NAME2(name,c20)), \
	bgq_params(NAME2(name,c21)), \
	bgq_params(NAME2(name,c22))

#define bgq_su3_spinor_decl(name) \
	bgq_su3_vdecl(NAME2(name,v0));     \
	bgq_su3_vdecl(NAME2(name,v1));     \
	bgq_su3_vdecl(NAME2(name,v2));     \
	bgq_su3_vdecl(NAME2(name,v3))

#define bgq_su3_spinor_vars(name) \
	bgq_su3_vvars(NAME2(name,v0)), \
	bgq_su3_vvars(NAME2(name,v1)), \
	bgq_su3_vvars(NAME2(name,v2)), \
	bgq_su3_vvars(NAME2(name,v3))

#define bgq_su3_weyl_vars(name) \
	bgq_su3_vvars(NAME2(name,v0)), \
	bgq_su3_vvars(NAME2(name,v1))

#define bgq_su3_spinor_params(name) \
	bgq_su3_vparams(NAME2(name,v0)), \
	bgq_su3_vparams(NAME2(name,v1)), \
	bgq_su3_vparams(NAME2(name,v2)), \
	bgq_su3_vparams(NAME2(name,v3))

#define bgq_su3_spinor_refparams(name) \
	bgq_su3_vrefparams(NAME2(name,v0)), \
	bgq_su3_vrefparams(NAME2(name,v1)), \
	bgq_su3_vrefparams(NAME2(name,v2)), \
	bgq_su3_vrefparams(NAME2(name,v3))

#define bgq_su3_weyl_decl(name) \
	bgq_su3_vdecl(NAME2(name,v0));	\
	bgq_su3_vdecl(NAME2(name,v1))

#define bgq_su3_weyl_decl_leftonly(name) \
	bgq_su3_vdecl_leftonly(NAME2(name,v0));	\
	bgq_su3_vdecl_leftonly(NAME2(name,v1))

#define bgq_su3_weyl_params(name) \
	bgq_su3_vparams(NAME2(name,v0)), \
	bgq_su3_vparams(NAME2(name,v1))

#define bgq_su3_vzero(dst) \
	bgq_zero(NAME2(dst,c0)); \
	bgq_zero(NAME2(dst,c1)); \
	bgq_zero(NAME2(dst,c2))

#define bgq_su3_spinor_zero(dst) \
	bgq_su3_vzero(NAME2(dst,v0));     \
	bgq_su3_vzero(NAME2(dst,v1));     \
	bgq_su3_vzero(NAME2(dst,v2));     \
	bgq_su3_vzero(NAME2(dst,v3))

#define bgq_su3_spinor_load NAME2(bgq_su3_spinor_load,PRECISION)
#if BGQ_QPX
#define bgq_su3_spinor_load_double(dst, addr) \
	do { \
		void *ptr = (addr); \
		bgq_qvlfduxa(NAME3(dst,v0,c0), ptr, 0); \
		bgq_qvlfduxa(NAME3(dst,v0,c1), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v0,c2), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v1,c0), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v1,c1), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v1,c2), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v2,c0), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v2,c1), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v2,c2), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v3,c0), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v3,c1), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v3,c2), ptr, 32); \
	} while (0)
#else
#define bgq_su3_spinor_load_double(dst, addr) \
	bgq_lda_double(NAME3(dst,v0,c0),   0, addr);          \
	bgq_lda_double(NAME3(dst,v0,c1),  32, addr);          \
	bgq_lda_double(NAME3(dst,v0,c2),  64, addr);          \
	bgq_lda_double(NAME3(dst,v1,c0),  96, addr);          \
	bgq_lda_double(NAME3(dst,v1,c1), 128, addr);          \
	bgq_lda_double(NAME3(dst,v1,c2), 160, addr);          \
	bgq_lda_double(NAME3(dst,v2,c0), 192, addr);          \
	bgq_lda_double(NAME3(dst,v2,c1), 224, addr);          \
	bgq_lda_double(NAME3(dst,v2,c2), 256, addr);          \
	bgq_lda_double(NAME3(dst,v3,c0), 288, addr);          \
	bgq_lda_double(NAME3(dst,v3,c1), 320, addr);          \
	bgq_lda_double(NAME3(dst,v3,c2), 352, addr)
#endif

#if BGQ_QPX
#define bgq_su3_spinor_load_float(dst, addr) \
	do { \
		void *ptr = (addr); \
		bgq_qvlfsuxa(NAME3(dst,v0,c0), ptr, 0); \
		bgq_qvlfsuxa(NAME3(dst,v0,c1), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v0,c2), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v1,c0), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v1,c1), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v1,c2), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v2,c0), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v2,c1), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v2,c2), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v3,c0), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v3,c1), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v3,c2), ptr, 16); \
	} while (0)
#else
#define bgq_su3_spinor_load_float(dst, addr) \
	bgq_lda_float(NAME3(dst,v0,c0),   0, addr);          \
	bgq_lda_float(NAME3(dst,v0,c1),  16, addr);          \
	bgq_lda_float(NAME3(dst,v0,c2),  32, addr);          \
	bgq_lda_float(NAME3(dst,v1,c0),  48, addr);          \
	bgq_lda_float(NAME3(dst,v1,c1),  64, addr);          \
	bgq_lda_float(NAME3(dst,v1,c2),  80, addr);          \
	bgq_lda_float(NAME3(dst,v2,c0),  96, addr);          \
	bgq_lda_float(NAME3(dst,v2,c1), 112, addr);          \
	bgq_lda_float(NAME3(dst,v2,c2), 128, addr);          \
	bgq_lda_float(NAME3(dst,v3,c0), 144, addr);          \
	bgq_lda_float(NAME3(dst,v3,c1), 160, addr);          \
	bgq_lda_float(NAME3(dst,v3,c2), 176, addr)
#endif

#define bgq_su3_weyl_load NAME2(bgq_su3_weyl_load,PRECISION)
#if BGQ_QPX
#define bgq_su3_weyl_load_double(dst, addr) \
	do { \
		void *ptr = (addr); \
		bgq_qvlfduxa(NAME3(dst,v0,c0), ptr, 0); \
		bgq_qvlfduxa(NAME3(dst,v0,c1), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v0,c2), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v1,c0), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v1,c1), ptr, 32); \
		bgq_qvlfduxa(NAME3(dst,v1,c2), ptr, 32); \
	} while (0)
#else
#define bgq_su3_weyl_load_double(dest, addr) \
	bgq_lda_double(NAME3(dest,v0,c0),   0, addr);        \
	bgq_lda_double(NAME3(dest,v0,c1),  32, addr);        \
	bgq_lda_double(NAME3(dest,v0,c2),  64, addr);        \
	bgq_lda_double(NAME3(dest,v1,c0),  96, addr);        \
	bgq_lda_double(NAME3(dest,v1,c1), 128, addr);        \
	bgq_lda_double(NAME3(dest,v1,c2), 160, addr)
#endif

#if BGQ_QPX
#define bgq_su3_weyl_load_float(dst, addr)       \
	do {                                         \
		void *ptr = (addr);                      \
		bgq_qvlfsuxa(NAME3(dst,v0,c0), ptr, 0);  \
		bgq_qvlfsuxa(NAME3(dst,v0,c1), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v0,c2), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v1,c0), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v1,c1), ptr, 16); \
		bgq_qvlfsuxa(NAME3(dst,v1,c2), ptr, 16); \
	} while (0)
#else
#define bgq_su3_weyl_load_float(dest, addr) \
	bgq_lda_float(NAME3(dest,v0,c0),   0, addr);        \
	bgq_lda_float(NAME3(dest,v0,c1),  16, addr);        \
	bgq_lda_float(NAME3(dest,v0,c2),  32, addr);        \
	bgq_lda_float(NAME3(dest,v1,c0),  48, addr);        \
	bgq_lda_float(NAME3(dest,v1,c1),  64, addr);        \
	bgq_lda_float(NAME3(dest,v1,c2),  80, addr)
#endif

// obsolete?
#define bgq_su3_weyl_left_move_double(dstaddr,srcaddr) \
	do { \
		bgq_vector4double_decl(tmp0); \
		bgq_vector4double_decl(tmp1); \
		bgq_vector4double_decl(tmp2); \
		bgq_lda_double(tmp0,  0, srcaddr); \
		bgq_lda_double(tmp1, 32, srcaddr); \
		bgq_lda_double(tmp2, 64, srcaddr); \
		bgq_sta_double(tmp0,  0, dstaddr); \
		bgq_sta_double(tmp1, 32, dstaddr); \
		bgq_sta_double(tmp2, 64, dstaddr); \
	} while (0)
	// 96 byte (3x 32B)

#define bgq_su3_weyl_load_left NAME2(bgq_su3_weyl_load_left,PRECISION)
#define bgq_su3_weyl_load_left_double(dest, addr) \
	bgq_ld2a_double(NAME3(dest,v0,c0),   0, addr);        \
	bgq_ld2a_double(NAME3(dest,v0,c1),  32, addr);        \
	bgq_ld2a_double(NAME3(dest,v0,c2),  64, addr);        \
	bgq_ld2a_double(NAME3(dest,v1,c0),  96, addr);        \
	bgq_ld2a_double(NAME3(dest,v1,c1), 128, addr);        \
	bgq_ld2a_double(NAME3(dest,v1,c2), 160, addr)

#define bgq_su3_weyl_load_left_float(dest, addr) \
	bgq_ld2a_float(NAME3(dest,v0,c0),   0, addr);        \
	bgq_ld2a_float(NAME3(dest,v0,c1),  16, addr);        \
	bgq_ld2a_float(NAME3(dest,v0,c2),  32, addr);        \
	bgq_ld2a_float(NAME3(dest,v1,c0),  48, addr);        \
	bgq_ld2a_float(NAME3(dest,v1,c1),  64, addr);        \
	bgq_ld2a_float(NAME3(dest,v1,c2),  80, addr)

#define bgq_su3_weyl_load_right NAME2(bgq_su3_weyl_load_right,PRECISION)
#define bgq_su3_weyl_load_right_double(dest, addr) \
	bgq_ld2a_double(NAME3(dest,v0,c0),   0+16, addr);        \
	bgq_ld2a_double(NAME3(dest,v0,c1),  32+16, addr);        \
	bgq_ld2a_double(NAME3(dest,v0,c2),  64+16, addr);        \
	bgq_ld2a_double(NAME3(dest,v1,c0),  96+16, addr);        \
	bgq_ld2a_double(NAME3(dest,v1,c1), 128+16, addr);        \
	bgq_ld2a_double(NAME3(dest,v1,c2), 160+16, addr)

#define bgq_su3_weyl_load_right_float(dest, addr) \
	bgq_ld2a_float(NAME3(dest,v0,c0),   0+8, addr);        \
	bgq_ld2a_float(NAME3(dest,v0,c1),  16+8, addr);        \
	bgq_ld2a_float(NAME3(dest,v0,c2),  32+8, addr);        \
	bgq_ld2a_float(NAME3(dest,v1,c0),  48+8, addr);        \
	bgq_ld2a_float(NAME3(dest,v1,c1),  64+8, addr);        \
	bgq_ld2a_float(NAME3(dest,v1,c2),  80+8, addr)

#define bgq_su3_matrix_valgen(dst) \
	do { \
		bgq_valgen(NAME2(dst,c00)); \
		bgq_valgen(NAME2(dst,c01)); \
		bgq_valgen(NAME2(dst,c02)); \
		bgq_valgen(NAME2(dst,c10)); \
		bgq_valgen(NAME2(dst,c11)); \
		bgq_valgen(NAME2(dst,c12)); \
		bgq_valgen(NAME2(dst,c20)); \
		bgq_valgen(NAME2(dst,c21)); \
		bgq_valgen(NAME2(dst,c22)); \
	} while (0)


#define bgq_su3_weyl_valgen(dst) \
	do { \
		bgq_valgen(NAME3(dst,v0,c0)); \
		bgq_valgen(NAME3(dst,v0,c1)); \
		bgq_valgen(NAME3(dst,v0,c2)); \
		bgq_valgen(NAME3(dst,v1,c0)); \
		bgq_valgen(NAME3(dst,v1,c1)); \
		bgq_valgen(NAME3(dst,v1,c2)); \
	} while (0)

#define bgq_su3_spinor_valgen(dst) \
	do { \
		bgq_valgen(NAME3(dst,v0,c0)); \
		bgq_valgen(NAME3(dst,v0,c1)); \
		bgq_valgen(NAME3(dst,v0,c2)); \
		bgq_valgen(NAME3(dst,v1,c0)); \
		bgq_valgen(NAME3(dst,v1,c1)); \
		bgq_valgen(NAME3(dst,v1,c2)); \
		bgq_valgen(NAME3(dst,v2,c0)); \
		bgq_valgen(NAME3(dst,v2,c1)); \
		bgq_valgen(NAME3(dst,v2,c2)); \
		bgq_valgen(NAME3(dst,v3,c0)); \
		bgq_valgen(NAME3(dst,v3,c1)); \
		bgq_valgen(NAME3(dst,v3,c2)); \
	} while (0)

#define bgq_su3_matrix_load NAME2(bgq_su3_matrix_load,PRECISION)
#if BGQ_QPX
#define bgq_su3_matrix_load_double(dst, addr) \
do {\
	void *ptr = (addr); \
	bgq_qvlfduxa(NAME2(dst,c00), ptr, 0); \
	bgq_qvlfduxa(NAME2(dst,c01), ptr, 32); \
	bgq_qvlfduxa(NAME2(dst,c02), ptr, 32); \
	bgq_qvlfduxa(NAME2(dst,c10), ptr, 32); \
	bgq_qvlfduxa(NAME2(dst,c11), ptr, 32); \
	bgq_qvlfduxa(NAME2(dst,c12), ptr, 32); \
	bgq_qvlfduxa(NAME2(dst,c20), ptr, 32); \
	bgq_qvlfduxa(NAME2(dst,c21), ptr, 32); \
	bgq_qvlfduxa(NAME2(dst,c22), ptr, 32); \
} while (0)
#else
#define bgq_su3_matrix_load_double(dest, addr) \
	bgq_lda_double(NAME2(dest,c00),   0, addr); \
	bgq_lda_double(NAME2(dest,c01),  32, addr); \
	bgq_lda_double(NAME2(dest,c02),  64, addr); \
	bgq_lda_double(NAME2(dest,c10),  96, addr); \
	bgq_lda_double(NAME2(dest,c11), 128, addr); \
	bgq_lda_double(NAME2(dest,c12), 160, addr); \
	bgq_lda_double(NAME2(dest,c20), 192, addr); \
	bgq_lda_double(NAME2(dest,c21), 224, addr); \
	bgq_lda_double(NAME2(dest,c22), 256, addr)
#endif

#define bgq_su3_matrix_load_float(dest, addr) \
	bgq_lda_float(NAME2(dest,c00),   0, addr); \
	bgq_lda_float(NAME2(dest,c01),  16, addr); \
	bgq_lda_float(NAME2(dest,c02),  32, addr); \
	bgq_lda_float(NAME2(dest,c10),  48, addr); \
	bgq_lda_float(NAME2(dest,c11),  64, addr); \
	bgq_lda_float(NAME2(dest,c12),  80, addr); \
	bgq_lda_float(NAME2(dest,c20),  96, addr); \
	bgq_lda_float(NAME2(dest,c21), 112, addr); \
	bgq_lda_float(NAME2(dest,c22), 128, addr)

#define bgq_su3_spinor_load_left NAME2(bgq_su3_spinor_load_left,PRECISION)
#define bgq_su3_spinor_load_left_double(dst,addr) \
	bgq_ld2a_double(NAME3(dst,v0,c0),   0, addr); \
	bgq_ld2a_double(NAME3(dst,v0,c1),  32, addr); \
	bgq_ld2a_double(NAME3(dst,v0,c2),  64, addr); \
	bgq_ld2a_double(NAME3(dst,v1,c0),  96, addr); \
	bgq_ld2a_double(NAME3(dst,v1,c1), 128, addr); \
	bgq_ld2a_double(NAME3(dst,v1,c2), 160, addr); \
	bgq_ld2a_double(NAME3(dst,v2,c0), 192, addr); \
	bgq_ld2a_double(NAME3(dst,v2,c1), 224, addr); \
	bgq_ld2a_double(NAME3(dst,v2,c2), 256, addr); \
	bgq_ld2a_double(NAME3(dst,v3,c0), 288, addr); \
	bgq_ld2a_double(NAME3(dst,v3,c1), 320, addr); \
	bgq_ld2a_double(NAME3(dst,v3,c2), 352, addr)

#define bgq_su3_spinor_load_left_float(dst,addr) \
	bgq_ld2a_float(NAME3(dst,v0,c0),   0, addr); \
	bgq_ld2a_float(NAME3(dst,v0,c1),  16, addr); \
	bgq_ld2a_float(NAME3(dst,v0,c2),  32, addr); \
	bgq_ld2a_float(NAME3(dst,v1,c0),  48, addr); \
	bgq_ld2a_float(NAME3(dst,v1,c1),  64, addr); \
	bgq_ld2a_float(NAME3(dst,v1,c2),  80, addr); \
	bgq_ld2a_float(NAME3(dst,v2,c0),  96, addr); \
	bgq_ld2a_float(NAME3(dst,v2,c1), 112, addr); \
	bgq_ld2a_float(NAME3(dst,v2,c2), 128, addr); \
	bgq_ld2a_float(NAME3(dst,v3,c0), 144, addr); \
	bgq_ld2a_float(NAME3(dst,v3,c1), 160, addr); \
	bgq_ld2a_float(NAME3(dst,v3,c2), 176, addr)

#define bgq_su3_spinor_load_right NAME2(bgq_su3_spinor_load_right,PRECISION)
#define bgq_su3_spinor_load_right_double(dst,addr) \
	bgq_ld2a_double(NAME3(dst,v0,c0),   0+16, addr); \
	bgq_ld2a_double(NAME3(dst,v0,c1),  32+16, addr); \
	bgq_ld2a_double(NAME3(dst,v0,c2),  64+16, addr); \
	bgq_ld2a_double(NAME3(dst,v1,c0),  96+16, addr); \
	bgq_ld2a_double(NAME3(dst,v1,c1), 128+16, addr); \
	bgq_ld2a_double(NAME3(dst,v1,c2), 160+16, addr); \
	bgq_ld2a_double(NAME3(dst,v2,c0), 192+16, addr); \
	bgq_ld2a_double(NAME3(dst,v2,c1), 224+16, addr); \
	bgq_ld2a_double(NAME3(dst,v2,c2), 256+16, addr); \
	bgq_ld2a_double(NAME3(dst,v3,c0), 288+16, addr); \
	bgq_ld2a_double(NAME3(dst,v3,c1), 320+16, addr); \
	bgq_ld2a_double(NAME3(dst,v3,c2), 352+16, addr)

#define bgq_su3_spinor_load_right_float(dst,addr) \
	bgq_ld2a_float(NAME3(dst,v0,c0),   0+8, addr); \
	bgq_ld2a_float(NAME3(dst,v0,c1),  16+8, addr); \
	bgq_ld2a_float(NAME3(dst,v0,c2),  32+8, addr); \
	bgq_ld2a_float(NAME3(dst,v1,c0),  48+8, addr); \
	bgq_ld2a_float(NAME3(dst,v1,c1),  64+8, addr); \
	bgq_ld2a_float(NAME3(dst,v1,c2),  80+8, addr); \
	bgq_ld2a_float(NAME3(dst,v2,c0),  96+8, addr); \
	bgq_ld2a_float(NAME3(dst,v2,c1), 112+8, addr); \
	bgq_ld2a_float(NAME3(dst,v2,c2), 128+8, addr); \
	bgq_ld2a_float(NAME3(dst,v3,c0), 144+8, addr); \
	bgq_ld2a_float(NAME3(dst,v3,c1), 160+8, addr); \
	bgq_ld2a_float(NAME3(dst,v3,c2), 176+8, addr)

#define bgq_su3_spinor_store NAME2(bgq_su3_spinor_store,PRECISION)
#if BGQ_QPX
#define bgq_su3_spinor_store_double(addr, spinor) \
	do {\
		void *ptr = (addr); \
		bgq_qvstfduxa(NAME3(spinor,v0,c0), ptr, 0); /*bgq_sta_double(NAME3(spinor,v0,c0),   0, addr);*/ /* qvstfdxa */ \
		bgq_qvstfduxa(NAME3(spinor,v0,c1), ptr, 32); \
		bgq_qvstfduxa(NAME3(spinor,v0,c2), ptr, 32); \
		bgq_qvstfduxa(NAME3(spinor,v1,c0), ptr, 32); \
		bgq_qvstfduxa(NAME3(spinor,v1,c1), ptr, 32); \
		bgq_qvstfduxa(NAME3(spinor,v1,c2), ptr, 32); \
		bgq_qvstfduxa(NAME3(spinor,v2,c0), ptr, 32); \
		bgq_qvstfduxa(NAME3(spinor,v2,c1), ptr, 32); \
		bgq_qvstfduxa(NAME3(spinor,v2,c2), ptr, 32); \
		bgq_qvstfduxa(NAME3(spinor,v3,c0), ptr, 32); \
		bgq_qvstfduxa(NAME3(spinor,v3,c1), ptr, 32); \
		bgq_qvstfduxa(NAME3(spinor,v3,c2), ptr, 32); \
	} while (0)
#else
#define bgq_su3_spinor_store_double(addr,src) \
	bgq_sta_double(NAME3(src,v0,c0),   0, addr);          \
	bgq_sta_double(NAME3(src,v0,c1),  32, addr);          \
	bgq_sta_double(NAME3(src,v0,c2),  64, addr);          \
	bgq_sta_double(NAME3(src,v1,c0),  96, addr);          \
	bgq_sta_double(NAME3(src,v1,c1), 128, addr);          \
	bgq_sta_double(NAME3(src,v1,c2), 160, addr);          \
	bgq_sta_double(NAME3(src,v2,c0), 192, addr);          \
	bgq_sta_double(NAME3(src,v2,c1), 224, addr);          \
	bgq_sta_double(NAME3(src,v2,c2), 256, addr);          \
	bgq_sta_double(NAME3(src,v3,c0), 288, addr);          \
	bgq_sta_double(NAME3(src,v3,c1), 320, addr);          \
	bgq_sta_double(NAME3(src,v3,c2), 352, addr)
#endif

#if BGQ_QPX
#define bgq_su3_spinor_store_float(addr,spinor) \
	do {\
		void *ptr = (addr); \
		bgq_qvstfsuxa(NAME3(spinor,v0,c0), ptr, 0);  \
		bgq_qvstfsuxa(NAME3(spinor,v0,c1), ptr, 16); \
		bgq_qvstfsuxa(NAME3(spinor,v0,c2), ptr, 16); \
		bgq_qvstfsuxa(NAME3(spinor,v1,c0), ptr, 16); \
		bgq_qvstfsuxa(NAME3(spinor,v1,c1), ptr, 16); \
		bgq_qvstfsuxa(NAME3(spinor,v1,c2), ptr, 16); \
		bgq_qvstfsuxa(NAME3(spinor,v2,c0), ptr, 16); \
		bgq_qvstfsuxa(NAME3(spinor,v2,c1), ptr, 16); \
		bgq_qvstfsuxa(NAME3(spinor,v2,c2), ptr, 16); \
		bgq_qvstfsuxa(NAME3(spinor,v3,c0), ptr, 16); \
		bgq_qvstfsuxa(NAME3(spinor,v3,c1), ptr, 16); \
		bgq_qvstfsuxa(NAME3(spinor,v3,c2), ptr, 16); \
	} while (0)
#else
#define bgq_su3_spinor_store_float(addr,src) \
	bgq_sta_float(NAME3(src,v0,c0),   0, addr);          \
	bgq_sta_float(NAME3(src,v0,c1),  16, addr);          \
	bgq_sta_float(NAME3(src,v0,c2),  32, addr);          \
	bgq_sta_float(NAME3(src,v1,c0),  48, addr);          \
	bgq_sta_float(NAME3(src,v1,c1),  64, addr);          \
	bgq_sta_float(NAME3(src,v1,c2),  80, addr);          \
	bgq_sta_float(NAME3(src,v2,c0),  96, addr);          \
	bgq_sta_float(NAME3(src,v2,c1), 112, addr);          \
	bgq_sta_float(NAME3(src,v2,c2), 128, addr);          \
	bgq_sta_float(NAME3(src,v3,c0), 144, addr);          \
	bgq_sta_float(NAME3(src,v3,c1), 160, addr);          \
	bgq_sta_float(NAME3(src,v3,c2), 176, addr)
#endif

#define bgq_su3_weyl_antioptaway(weyl) \
	do { \
		bgq_antioptaway(NAME3(weyl,v0,c0)); \
		bgq_antioptaway(NAME3(weyl,v0,c1)); \
		bgq_antioptaway(NAME3(weyl,v0,c2)); \
		bgq_antioptaway(NAME3(weyl,v1,c0)); \
		bgq_antioptaway(NAME3(weyl,v1,c1)); \
		bgq_antioptaway(NAME3(weyl,v1,c2)); \
	} while (0)

#define bgq_su3_weyl_store NAME2(bgq_su3_weyl_store,PRECISION)
#if BGQ_QPX
#define bgq_su3_weyl_store_double(addr,weyl) \
	do { \
		void *ptr = (addr); \
		bgq_qvstfduxa(NAME3(weyl,v0,c0), ptr, 0); \
		bgq_qvstfduxa(NAME3(weyl,v0,c1), ptr, 32); \
		bgq_qvstfduxa(NAME3(weyl,v0,c2), ptr, 32); \
		bgq_qvstfduxa(NAME3(weyl,v1,c0), ptr, 32); \
		bgq_qvstfduxa(NAME3(weyl,v1,c1), ptr, 32); \
		bgq_qvstfduxa(NAME3(weyl,v1,c2), ptr, 32); \
	} while (0)
#else
#define bgq_su3_weyl_store_double(addr,src) \
	bgq_sta_double(NAME3(src,v0,c0),   0, addr);          \
	bgq_sta_double(NAME3(src,v0,c1),  32, addr);          \
	bgq_sta_double(NAME3(src,v0,c2),  64, addr);          \
	bgq_sta_double(NAME3(src,v1,c0),  96, addr);          \
	bgq_sta_double(NAME3(src,v1,c1), 128, addr);          \
	bgq_sta_double(NAME3(src,v1,c2), 160, addr)
	/* 196 bytes = 3 L1 lines */
#endif


#if BGQ_QPX
#define bgq_su3_weyl_store_float(addr,weyl) \
	do { \
		void *ptr = (addr); \
		bgq_qvstfsuxa(NAME3(weyl,v0,c0), ptr, 0); \
		bgq_qvstfsuxa(NAME3(weyl,v0,c1), ptr, 16); \
		bgq_qvstfsuxa(NAME3(weyl,v0,c2), ptr, 16); \
		bgq_qvstfsuxa(NAME3(weyl,v1,c0), ptr, 16); \
		bgq_qvstfsuxa(NAME3(weyl,v1,c1), ptr, 16); \
		bgq_qvstfsuxa(NAME3(weyl,v1,c2), ptr, 16); \
	} while (0)
#else
#define bgq_su3_weyl_store_float(addr,src) \
	bgq_sta_float(NAME3(src,v0,c0),   0, addr);          \
	bgq_sta_float(NAME3(src,v0,c1),  16, addr);          \
	bgq_sta_float(NAME3(src,v0,c2),  32, addr);          \
	bgq_sta_float(NAME3(src,v1,c0),  48, addr);          \
	bgq_sta_float(NAME3(src,v1,c1),  64, addr);          \
	bgq_sta_float(NAME3(src,v1,c2),  80, addr)
#endif

#define bgq_su3_weyl_left_store_double(addr,src) \
	do { \
		void *ptr = (addr); \
		bgq_st2a_double(NAME3(src,v0,c0), 0, ptr); \
		bgq_qvstfcduxa(NAME3(src,v0,c1), ptr, 16); \
		bgq_qvstfcduxa(NAME3(src,v0,c2), ptr, 16); \
		bgq_qvstfcduxa(NAME3(src,v1,c0), ptr, 16); \
		bgq_qvstfcduxa(NAME3(src,v1,c1), ptr, 16); \
		bgq_qvstfcduxa(NAME3(src,v1,c2), ptr, 16); \
	} while (0)
	// 2*3*2*8 = 96 byte

#define bgq_su3_weyl_right_store_double(addr,src) \
	do { \
		void *ptr = (addr); \
		bgq_vector4double_decl(right0); \
		bgq_vector4double_decl(right1); \
		bgq_vector4double_decl(right2); \
		bgq_rmerge(right0, NAME3(src,v0,c0), NAME3(src,v0,c1)); \
		bgq_rmerge(right1, NAME3(src,v0,c2), NAME3(src,v1,c0)); \
		bgq_rmerge(right2, NAME3(src,v1,c1), NAME3(src,v1,c2)); \
		bgq_sta_double(right0, 0, ptr); \
		bgq_qvstfduxa(right1, ptr, 32); \
		bgq_qvstfduxa(right2, ptr, 32); \
	} while (0)
// 96 byte

#define bgq_su3_spinor_merge(dst,a,b)         \
	bgq_su3_vmerge(dst##_v0, a##_v0, b##_v0); \
	bgq_su3_vmerge(dst##_v1, a##_v1, b##_v1); \
	bgq_su3_vmerge(dst##_v2, a##_v2, b##_v2); \
	bgq_su3_vmerge(dst##_v3, a##_v3, b##_v3)

#define bgq_su3_weyl_merge2(dst,a,b)         \
	bgq_su3_vmerge2(NAME2(dst,v0), NAME2(a,v0), NAME2(b,v0)); \
	bgq_su3_vmerge2(NAME2(dst,v1), NAME2(a,v1), NAME2(b,v1))

#define bgq_su3_weyl_lmerge(dst,a,b)         \
	bgq_su3_vlmerge(NAME2(dst,v0), NAME2(a,v0), NAME2(b,v0)); \
	bgq_su3_vlmerge(NAME2(dst,v1), NAME2(a,v1), NAME2(b,v1))

#define bgq_su3_weyl_rmerge(dst,a,b)         \
	bgq_su3_vrmerge(NAME2(dst,v0), NAME2(a,v0), NAME2(b,v0)); \
	bgq_su3_vrmerge(NAME2(dst,v1), NAME2(a,v1), NAME2(b,v1))

#define bgq_su3_vmov(dst,src)    \
	bgq_mov(NAME2(dst,c0), NAME2(src,c0)); \
	bgq_mov(NAME2(dst,c1), NAME2(src,c1)); \
	bgq_mov(NAME2(dst,c2), NAME2(src,c2))

#define bgq_su3_vmerge2(dst,a,b)           \
	bgq_merge2(NAME2(dst,c0), NAME2(a,c0), NAME2(b,c0)); \
	bgq_merge2(NAME2(dst,c1), NAME2(a,c1), NAME2(b,c1)); \
	bgq_merge2(NAME2(dst,c2), NAME2(a,c2), NAME2(b,c2))

#define bgq_su3_vlmerge(dst,a,b)           \
	bgq_lmerge(NAME2(dst,c0), NAME2(a,c0), NAME2(b,c0)); \
	bgq_lmerge(NAME2(dst,c1), NAME2(a,c1), NAME2(b,c1)); \
	bgq_lmerge(NAME2(dst,c2), NAME2(a,c2), NAME2(b,c2))

#define bgq_su3_vrmerge(dst,a,b)           \
	bgq_rmerge(NAME2(dst,c0), NAME2(a,c0), NAME2(b,c0)); \
	bgq_rmerge(NAME2(dst,c1), NAME2(a,c1), NAME2(b,c1)); \
	bgq_rmerge(NAME2(dst,c2), NAME2(a,c2), NAME2(b,c2))

#define bgq_su3_mmerge(dst,a,b)               \
	dst##_c00 = cvec_merge(a##_c00, b##_c00); \
	dst##_c01 = cvec_merge(a##_c01, b##_c01); \
	dst##_c02 = cvec_merge(a##_c02, b##_c02); \
	dst##_c10 = cvec_merge(a##_c10, b##_c10); \
	dst##_c11 = cvec_merge(a##_c11, b##_c11); \
	dst##_c12 = cvec_merge(a##_c12, b##_c12); \
	dst##_c20 = cvec_merge(a##_c20, b##_c20); \
	dst##_c21 = cvec_merge(a##_c21, b##_c21); \
	dst##_c22 = cvec_merge(a##_c22, b##_c22); \

#define bgq_su3_vadd(dst,v1,v2)          \
	bgq_add(NAME2(dst,c0), NAME2(v1,c0), NAME2(v2,c0)); \
	bgq_add(NAME2(dst,c1), NAME2(v1,c1), NAME2(v2,c1)); \
	bgq_add(NAME2(dst,c2), NAME2(v1,c2), NAME2(v2,c2))

#define bgq_su3_vsub(dst,v1,v2)         \
	bgq_sub(NAME2(dst,c0), NAME2(v1,c0), NAME2(v2,c0)); \
	bgq_sub(NAME2(dst,c1), NAME2(v1,c1), NAME2(v2,c1)); \
	bgq_sub(NAME2(dst,c2), NAME2(v1,c2), NAME2(v2,c2))

#define bgq_su3_vpiadd(dst,v1,v2)         \
	bgq_iadd(NAME2(dst,c0), NAME2(v1,c0), NAME2(v2,c0)); \
	bgq_iadd(NAME2(dst,c1), NAME2(v1,c1), NAME2(v2,c1)); \
	bgq_iadd(NAME2(dst,c2), NAME2(v1,c2), NAME2(v2,c2))

#define bgq_su3_vpisub(dst,v1,v2)         \
	bgq_isub(NAME2(dst,c0), NAME2(v1,c0), NAME2(v2,c0)); \
	bgq_isub(NAME2(dst,c1), NAME2(v1,c1), NAME2(v2,c1)); \
	bgq_isub(NAME2(dst,c2), NAME2(v1,c2), NAME2(v2,c2))

#define bgq_su3_vmul(dst,r,v) \
	do { \
		bgq_mul(NAME2(dst,c0), r, NAME2(v,c0)); \
		bgq_mul(NAME2(dst,c1), r, NAME2(v,c1)); \
		bgq_mul(NAME2(dst,c2), r, NAME2(v,c2)); \
	} while (0)

#define bgq_su3_cvmul(dst,c,v)              \
	bgq_cmul(NAME2(dst,c0), c, NAME2(v,c0)); \
	bgq_cmul(NAME2(dst,c1), c, NAME2(v,c1)); \
	bgq_cmul(NAME2(dst,c2), c, NAME2(v,c2))


#define bgq_su3_cvmadd(dst,c,v,vadd)              \
	do { \
		bgq_cmadd(NAME2(dst,c0), c, NAME2(v,c0), NAME2(vadd,c0)); \
		bgq_cmadd(NAME2(dst,c1), c, NAME2(v,c1), NAME2(vadd,c0)); \
		bgq_cmadd(NAME2(dst,c2), c, NAME2(v,c2), NAME2(vadd,c0)); \
	} while (0)

#define bgq_su3_cjgvmul(dst,c,v)              \
	bgq_cjgmul(NAME2(dst,c0), c, NAME2(v,c0)); \
	bgq_cjgmul(NAME2(dst,c1), c, NAME2(v,c1)); \
	bgq_cjgmul(NAME2(dst,c2), c, NAME2(v,c2))


#define bgq_su3_weyl_cmul(dst,c,weyl) \
	bgq_su3_cvmul(NAME2(dst,v0), c, NAME2(weyl,v0)); \
	bgq_su3_cvmul(NAME2(dst,v1), c, NAME2(weyl,v1))

#define bgq_su3_weyl_cjgmul(dst,c,weyl) \
	bgq_su3_cjgvmul(NAME2(dst,v0), c, NAME2(weyl,v0)); \
	bgq_su3_cjgvmul(NAME2(dst,v1), c, NAME2(weyl,v1))


#define bgq_su3_spinor_mul(dst,r,spinor) \
	do { \
		bgq_su3_vmul(NAME2(dst,v0),r,NAME2(spinor,v0)); \
		bgq_su3_vmul(NAME2(dst,v1),r,NAME2(spinor,v1)); \
		bgq_su3_vmul(NAME2(dst,v2),r,NAME2(spinor,v2)); \
		bgq_su3_vmul(NAME2(dst,v3),r,NAME2(spinor,v3)); \
	} while (0)

#define bgq_su3_spinor_add(dst,lhs,rhs) \
	do { \
		bgq_su3_vadd(NAME2(dst,v0),NAME2(lhs,v0),NAME2(rhs,v0)); \
		bgq_su3_vadd(NAME2(dst,v1),NAME2(lhs,v1),NAME2(rhs,v1)); \
		bgq_su3_vadd(NAME2(dst,v2),NAME2(lhs,v2),NAME2(rhs,v2)); \
		bgq_su3_vadd(NAME2(dst,v3),NAME2(lhs,v3),NAME2(rhs,v3)); \
	} while (0)


#define bgq_su3_mvmul(dst,m,v)                                                      \
	{                                                                               \
		bgq_su3_vdecl(MAKENAME4(mvmul,dst,m,v));                           \
		bgq_cmul (MAKENAME5(mvmul,dst,m,v,c0), NAME2(v,c0), NAME2(m,c00));         \
		bgq_cmadd(MAKENAME5(mvmul,dst,m,v,c0), NAME2(v,c1), NAME2(m,c01), MAKENAME5(mvmul,dst,m,v,c0)); \
		bgq_cmadd(MAKENAME5(mvmul,dst,m,v,c0), NAME2(v,c2), NAME2(m,c02), MAKENAME5(mvmul,dst,m,v,c0)); \
		bgq_cmul (MAKENAME5(mvmul,dst,m,v,c1), NAME2(v,c0), NAME2(m,c10));         \
		bgq_cmadd(MAKENAME5(mvmul,dst,m,v,c1), NAME2(v,c1), NAME2(m,c11), MAKENAME5(mvmul,dst,m,v,c1)); \
		bgq_cmadd(MAKENAME5(mvmul,dst,m,v,c1), NAME2(v,c2), NAME2(m,c12), MAKENAME5(mvmul,dst,m,v,c1)); \
		bgq_cmul (MAKENAME5(mvmul,dst,m,v,c2), NAME2(v,c0), NAME2(m,c20));         \
		bgq_cmadd(MAKENAME5(mvmul,dst,m,v,c2), NAME2(v,c1), NAME2(m,c21), MAKENAME5(mvmul,dst,m,v,c2)); \
		bgq_cmadd(MAKENAME5(mvmul,dst,m,v,c2), NAME2(v,c2), NAME2(m,c22), MAKENAME5(mvmul,dst,m,v,c2)); \
		bgq_su3_vmov(dst, MAKENAME4(mvmul,dst,m,v));                       \
	}

#define bgq_su3_mvinvmul(dst,m,v)                    \
	{                                                \
		bgq_su3_vdecl(MAKENAME4(mvmul,dst,m,v));                          \
		bgq_ccmul (MAKENAME5(mvmul,dst,m,v,c0), NAME2(v,c0), NAME2(m,c00));         \
		bgq_ccmadd(MAKENAME5(mvmul,dst,m,v,c0), NAME2(v,c1), NAME2(m,c10), MAKENAME5(mvmul,dst,m,v,c0)); \
		bgq_ccmadd(MAKENAME5(mvmul,dst,m,v,c0), NAME2(v,c2), NAME2(m,c20), MAKENAME5(mvmul,dst,m,v,c0)); \
		bgq_ccmul (MAKENAME5(mvmul,dst,m,v,c1), NAME2(v,c0), NAME2(m,c01));         \
		bgq_ccmadd(MAKENAME5(mvmul,dst,m,v,c1), NAME2(v,c1), NAME2(m,c11), MAKENAME5(mvmul,dst,m,v,c1)); \
		bgq_ccmadd(MAKENAME5(mvmul,dst,m,v,c1), NAME2(v,c2), NAME2(m,c21), MAKENAME5(mvmul,dst,m,v,c1)); \
		bgq_ccmul (MAKENAME5(mvmul,dst,m,v,c2), NAME2(v,c0), NAME2(m,c02));         \
		bgq_ccmadd(MAKENAME5(mvmul,dst,m,v,c2), NAME2(v,c1), NAME2(m,c12), MAKENAME5(mvmul,dst,m,v,c2)); \
		bgq_ccmadd(MAKENAME5(mvmul,dst,m,v,c2), NAME2(v,c2), NAME2(m,c22), MAKENAME5(mvmul,dst,m,v,c2)); \
		bgq_su3_vmov(dst, MAKENAME4(mvmul,dst,m,v));                      \
	}

#define bgq_su3_weyl_mvmul(dst,m,weyl)              \
	bgq_su3_mvmul(NAME2(dst,v0), m, NAME2(weyl,v0)); \
	bgq_su3_mvmul(NAME2(dst,v1), m, NAME2(weyl,v1))

#define bgq_su3_weyl_mvinvmul(dst,m,weyl) \
	bgq_su3_mvinvmul(NAME2(dst,v0), m, NAME2(weyl,v0)); \
	bgq_su3_mvinvmul(NAME2(dst,v1), m, NAME2(weyl,v1))


#define bgq_su3_spinor_mov(dst,src)  \
	bgq_su3_vmov(NAME2(dst,v0),NAME2(src,v0)); \
	bgq_su3_vmov(NAME2(dst,v1),NAME2(src,v1)); \
	bgq_su3_vmov(NAME2(dst,v2),NAME2(src,v2)); \
	bgq_su3_vmov(NAME2(dst,v3),NAME2(src,v3))

#define bgq_su3_weyl_mov(dst,src)  \
	bgq_su3_vmov(NAME2(dst,v0),NAME2(src,v0)); \
	bgq_su3_vmov(NAME2(dst,v1),NAME2(src,v1))


#define bgq_su3_reduce_weyl_tup(weyl, spinor) \
	bgq_su3_vadd(NAME2(weyl,v0), NAME2(spinor,v0), NAME2(spinor,v2)); \
	bgq_su3_vadd(NAME2(weyl,v1), NAME2(spinor,v1), NAME2(spinor,v3))

#define bgq_su3_reduce_weyl_tdown(weyl, spinor) \
	bgq_su3_vsub(NAME2(weyl,v0), NAME2(spinor,v0), NAME2(spinor,v2)); \
	bgq_su3_vsub(NAME2(weyl,v1), NAME2(spinor,v1), NAME2(spinor,v3))

#define bgq_su3_reduce_weyl_xup(weyl, spinor) \
	bgq_su3_vpiadd(NAME2(weyl,v0), NAME2(spinor,v0), NAME2(spinor,v3)); \
	bgq_su3_vpiadd(NAME2(weyl,v1), NAME2(spinor,v1), NAME2(spinor,v2))

#define bgq_su3_reduce_weyl_xdown(weyl, spinor) \
	bgq_su3_vpisub(NAME2(weyl,v0), NAME2(spinor,v0), NAME2(spinor,v3)); \
	bgq_su3_vpisub(NAME2(weyl,v1), NAME2(spinor,v1), NAME2(spinor,v2))

#define bgq_su3_reduce_weyl_yup(weyl, spinor) \
	bgq_su3_vadd(NAME2(weyl,v0), NAME2(spinor,v0), NAME2(spinor,v3)); \
	bgq_su3_vsub(NAME2(weyl,v1), NAME2(spinor,v1), NAME2(spinor,v2))

#define bgq_su3_reduce_weyl_ydown(weyl, spinor) \
	bgq_su3_vsub(NAME2(weyl,v0), NAME2(spinor,v0), NAME2(spinor,v3)); \
	bgq_su3_vadd(NAME2(weyl,v1), NAME2(spinor,v1), NAME2(spinor,v2))

#define bgq_su3_reduce_weyl_zup(weyl, spinor) \
	bgq_su3_vpiadd(NAME2(weyl,v0), NAME2(spinor,v0), NAME2(spinor,v2)); \
	bgq_su3_vpisub(NAME2(weyl,v1), NAME2(spinor,v1), NAME2(spinor,v3))

#define bgq_su3_reduce_weyl_zdown(weyl, spinor) \
	bgq_su3_vpisub(NAME2(weyl,v0), NAME2(spinor,v0), NAME2(spinor,v2)); \
	bgq_su3_vpiadd(NAME2(weyl,v1), NAME2(spinor,v1), NAME2(spinor,v3))


#define bgq_su3_expand_weyl_tup(result, weyl) \
	bgq_su3_vmov(NAME2(result,v0), NAME2(weyl,v0)); \
	bgq_su3_vmov(NAME2(result,v1), NAME2(weyl,v1)); \
	bgq_su3_vmov(NAME2(result,v2), NAME2(weyl,v0)); \
	bgq_su3_vmov(NAME2(result,v3), NAME2(weyl,v1))

#define bgq_su3_accum_weyl_tup(result, weyl) \
	bgq_su3_vadd(NAME2(result,v0), NAME2(result,v0), NAME2(weyl,v0)); \
	bgq_su3_vadd(NAME2(result,v1), NAME2(result,v1), NAME2(weyl,v1)); \
	bgq_su3_vadd(NAME2(result,v2), NAME2(result,v2), NAME2(weyl,v0)); \
	bgq_su3_vadd(NAME2(result,v3), NAME2(result,v3), NAME2(weyl,v1))

#define bgq_su3_accum_weyl_tdown(result, weyl) \
	bgq_su3_vadd(NAME2(result,v0), NAME2(result,v0), NAME2(weyl,v0)); \
	bgq_su3_vadd(NAME2(result,v1), NAME2(result,v1), NAME2(weyl,v1)); \
	bgq_su3_vsub(NAME2(result,v2), NAME2(result,v2), NAME2(weyl,v0)); \
	bgq_su3_vsub(NAME2(result,v3), NAME2(result,v3), NAME2(weyl,v1))

#define bgq_su3_accum_weyl_xup(result, weyl) \
	bgq_su3_vadd(NAME2(result,v0), NAME2(result,v0), NAME2(weyl,v0)); \
	bgq_su3_vadd(NAME2(result,v1), NAME2(result,v1), NAME2(weyl,v1)); \
	bgq_su3_vpisub(NAME2(result,v2), NAME2(result,v2), NAME2(weyl,v1)); \
	bgq_su3_vpisub(NAME2(result,v3), NAME2(result,v3), NAME2(weyl,v0))

#define bgq_su3_accum_weyl_xdown(result, weyl) \
	bgq_su3_vadd(NAME2(result,v0), NAME2(result,v0), NAME2(weyl,v0)); \
	bgq_su3_vadd(NAME2(result,v1), NAME2(result,v1), NAME2(weyl,v1)); \
	bgq_su3_vpiadd(NAME2(result,v2), NAME2(result,v2), NAME2(weyl,v1)); \
	bgq_su3_vpiadd(NAME2(result,v3), NAME2(result,v3), NAME2(weyl,v0))

#define bgq_su3_accum_weyl_yup(result, weyl) \
	bgq_su3_vadd(NAME2(result,v0), NAME2(result,v0), NAME2(weyl,v0)); \
	bgq_su3_vadd(NAME2(result,v1), NAME2(result,v1), NAME2(weyl,v1)); \
	bgq_su3_vsub(NAME2(result,v2), NAME2(result,v2), NAME2(weyl,v1)); \
	bgq_su3_vadd(NAME2(result,v3), NAME2(result,v3), NAME2(weyl,v0))

#define bgq_su3_accum_weyl_ydown(result, weyl) \
	bgq_su3_vadd(NAME2(result,v0), NAME2(result,v0), NAME2(weyl,v0)); \
	bgq_su3_vadd(NAME2(result,v1), NAME2(result,v1), NAME2(weyl,v1)); \
	bgq_su3_vadd(NAME2(result,v2), NAME2(result,v2), NAME2(weyl,v1)); \
	bgq_su3_vsub(NAME2(result,v3), NAME2(result,v3), NAME2(weyl,v0))

#define bgq_su3_accum_weyl_zup(result, weyl) \
	bgq_su3_vadd(NAME2(result,v0), NAME2(result,v0), NAME2(weyl,v0)); \
	bgq_su3_vadd(NAME2(result,v1), NAME2(result,v1), NAME2(weyl,v1)); \
	bgq_su3_vpisub(NAME2(result,v2), NAME2(result,v2), NAME2(weyl,v0)); \
	bgq_su3_vpiadd(NAME2(result,v3), NAME2(result,v3), NAME2(weyl,v1))

#define bgq_su3_accum_weyl_zdown(result, weyl) \
	bgq_su3_vadd(NAME2(result,v0), NAME2(result,v0), NAME2(weyl,v0)); \
	bgq_su3_vadd(NAME2(result,v1), NAME2(result,v1), NAME2(weyl,v1)); \
	bgq_su3_vpiadd(NAME2(result,v2), NAME2(result,v2), NAME2(weyl,v0)); \
	bgq_su3_vpisub(NAME2(result,v3), NAME2(result,v3), NAME2(weyl,v1))






#if !BGQ_QPX

// No semantic effects
#define bgq_prefetch(addr)
#define bgq_prefetchforwrite(addr)
#define bgq_prefetch_forward(addr)
#define bgq_prefetch_backward(addr)
#define bgq_flush(addr)

#define bgq_l1_zero(addr) \
	/*memset((addr),0,64)*/

#else

#if defined(XLC)
	#define bgq_prefetch(addr) \
		__dcbt(addr)
// __prefetch_by_load(addr) generates an lbz instruction, which is 'blocking'
// #include <spi/mu/Util.h>
// muspi_dcbt(vaddr,inc)
// muspi_dcbt_l2(vaddr,inc)
	#define bgq_prefetchforwrite(addr) \
		__dcbtst(addr)
	#define bgq_prefetch_forward(addr) \
		__prefetch_by_stream(1/*forward*/,(addr))
	#define bgq_prefetch_backward(addr) \
		__prefetch_by_stream(3/*backward*/,(addr))
	#define bgq_l1_zero(addr) \
		__dcbz(addr) /* sets 128 bytes (L2 cache line size) to zero */
//#define muspi_dcbz(vaddr,inc)
	#define bgq_flush(addr) \
		__dcbf(addr)
#elif defined(__GNUC__)
	#define bgq_prefetch(addr) \
		__builtin_prefetch((addr),0/*read*/)
	#define bgq_prefetchforwrite(addr) \
		__builtin_prefetch((addr),1/*write*/)
	#define bgq_prefetch_forward(addr) \
		bgq_prefetch(addr)
	#define bgq_prefetch_backward(addr) \
		bgq_prefetch(addr)
	#define bgq_l1_zero(addr)
	#define bgq_flush(addr)
#else
	#define bgq_prefetch(addr)
	#define bgq_prefetchforwrite(addr)
	#define bgq_prefetch_forward(addr)
	#define bgq_prefetch_backward(addr)
	#define bgq_l1_zero(addr)
	#define bgq_flush(addr)
#endif

#endif


#define bgq_ptrnext_prefetch(addr) \
	bgq_prefetch((void**)(addr) + 1)



#define bgq_su3_spinor_prefetch NAME2(bgq_su3_spinor_prefetch,PRECISION)
#if BGQ_QPX
#define bgq_su3_spinor_prefetch_double(addr)     \
	do { \
		void *ptr = (addr); \
		asm ( \
			"dcbt       0,%[ptr]  \n" \
			"dcbt  %[c64],%[ptr]  \n" \
			"dcbt %[c128],%[ptr]  \n" \
			"dcbt %[c192],%[ptr]  \n" \
			"dcbt %[c256],%[ptr]  \n" \
			"dcbt %[c320],%[ptr]  \n" \
			: : \
	 	 	  [ptr] "+r" (ptr), \
			  [c64]  "b" (64),  \
			  [c128] "b" (128), \
			  [c192] "b" (192), \
			  [c256] "b" (256), \
			  [c320] "b" (320)  \
		); \
	} while (0)
#else
#define bgq_su3_spinor_prefetch_double(addr)     \
	bgq_prefetch((char*)(addr) +   0);    \
	bgq_prefetch((char*)(addr) +  64);    \
	bgq_prefetch((char*)(addr) + 128);    \
	bgq_prefetch((char*)(addr) + 192);    \
	bgq_prefetch((char*)(addr) + 256);    \
	bgq_prefetch((char*)(addr) + 320)
#endif

#if BGQ_QPX
#define bgq_su3_spinor_prefetch_float(addr)     \
	do { \
		void *ptr = (addr); \
		asm ( \
			"dcbt       0,%[ptr]  \n" \
			"dcbt  %[c64],%[ptr]  \n" \
			"dcbt %[c128],%[ptr]  \n" \
			: : \
			  [ptr]  "r" (ptr), \
			  [c64]  "b" (64), \
			  [c128] "b" (128) \
		); \
	} while (0)
#else
#define bgq_su3_spinor_prefetch_float(addr)     \
	bgq_prefetch((char*)(addr) +   0);    \
	bgq_prefetch((char*)(addr) +  64);    \
	bgq_prefetch((char*)(addr) + 128)
#endif


#define bgq_su3_weyl_prefetch NAME2(bgq_su3_weyl_prefetch,PRECISION)
#if BGQ_QPX
#define bgq_su3_weyl_prefetch_double(addr)      \
	do { \
		void *ptr = (addr); \
		asm ( \
			"dcbt      0 ,%[ptr]  \n" \
			"dcbt  %[c64],%[ptr]  \n" \
			"dcbt %[c128],%[ptr]  \n" \
			: : \
			 [ptr] "r" (ptr), \
			 [c64] "b" (64), \
			 [c128] "b" (128) \
		); \
	} while (0)
#else
#define bgq_su3_weyl_prefetch_double(addr)      \
	bgq_prefetch((char*)(addr) +   0);    \
	bgq_prefetch((char*)(addr) +  64);    \
	bgq_prefetch((char*)(addr) + 128)
#endif

#define bgq_su3_weylnext_prefetch NAME2(bgq_su3_weylnext_prefetch,PRECISION)
#if BGQ_QPX
#define bgq_su3_weylnext_prefetch_double(addr)     \
	do {                             \
		void *ptr = (addr);          \
		asm (                        \
			"dcbt   %[c0],%[ptr]  \n" \
			"dcbt  %[c64],%[ptr]  \n" \
			"dcbt %[c128],%[ptr]  \n" \
			: : [ptr] "r" (ptr),       \
			    [c0] "b" (192+0),     \
			   [c64] "b" (192+64),    \
			  [c128] "b" (192+128)    \
		);                            \
	} while (0)
#else
#define bgq_su3_weylnext_prefetch_double(addr)      \
	bgq_prefetch((char*)(addr) + 192 +   0);    \
	bgq_prefetch((char*)(addr) + 192 +  64);    \
	bgq_prefetch((char*)(addr) + 192 + 128)
#endif

#if BGQ_QPX
#define bgq_su3_weylnext_prefetch_float(addr)     \
	do {                             \
		void *ptr = (addr);          \
		asm (                        \
			"dcbt   %[c0],%[ptr]  \n" \
			"dcbt  %[c64],%[ptr]  \n" \
			: : [ptr] "r" (ptr),       \
			    [c0] "b" (96+0),      \
			   [c64] "b" (96+64)      \
		);                            \
	} while (0)
#else
#define bgq_su3_weylnext_prefetch_float(addr)      \
	bgq_prefetch((char*)(addr) + 192 +   0);    \
	bgq_prefetch((char*)(addr) + 192 +  64);
#endif



#define bgq_su3_spinornext_prefetch NAME2(bgq_su3_spinornext_prefetch,PRECISION)
#if BGQ_QPX
#define bgq_su3_spinornext_prefetch_double(addr)     \
	do {                             \
		void *ptr = (addr);          \
		asm (                        \
			"dcbt   %[c0],%[ptr]  \n" \
			"dcbt  %[c64],%[ptr]  \n" \
			"dcbt %[c128],%[ptr]  \n" \
			"dcbt %[c192],%[ptr]  \n" \
			"dcbt %[c256],%[ptr]  \n" \
			"dcbt %[c320],%[ptr]  \n" \
			: : [ptr] "r" (ptr),       \
			    [c0] "b" (384+0),     \
			   [c64] "b" (384+64),    \
			  [c128] "b" (384+128),   \
			  [c192] "b" (384+192),   \
			  [c256] "b" (384+256),   \
			  [c320] "b" (384+320)   \
		);                            \
	} while (0)
#else
#define bgq_su3_spinornext_prefetch_double(addr)      \
	bgq_prefetch((char*)(addr) + 384 +   0);    \
	bgq_prefetch((char*)(addr) + 384 +  64);    \
	bgq_prefetch((char*)(addr) + 384 + 128);    \
	bgq_prefetch((char*)(addr) + 384 + 192);    \
	bgq_prefetch((char*)(addr) + 384 + 256);    \
	bgq_prefetch((char*)(addr) + 384 + 320)
#endif

#if BGQ_QPX
#define bgq_su3_spinornext_prefetch_float(addr)     \
	do {                             \
		void *ptr = (addr);          \
		asm (                        \
			"dcbt   %[c0],%[ptr]  \n" \
			"dcbt  %[c64],%[ptr]  \n" \
			"dcbt %[c128],%[ptr]  \n" \
			: : [ptr] "r" (ptr),       \
				[c0] "b" (192+0),     \
			   [c64] "b" (192+64),    \
			  [c128] "b" (192+128)   \
		);                            \
	} while (0)
#else
#define bgq_su3_spinornext_prefetch_float(addr)      \
	bgq_prefetch((char*)(addr) + 192 +   0);    \
	bgq_prefetch((char*)(addr) + 192 +  64);    \
	bgq_prefetch((char*)(addr) + 192 + 128)
#endif


#define bgq_su3_weylnextnext_prefetch NAME2(bgq_su3_weylnextnext_prefetch,PRECISION)
#if BGQ_QPX
#define bgq_su3_weylnextnext_prefetch_double(addr)     \
	do {                             \
		void *ptr = (addr);          \
		asm (                        \
			"dcbt   %[c0],%[ptr]  \n" \
			"dcbt  %[c64],%[ptr]  \n" \
			"dcbt %[c128],%[ptr]  \n" \
			: : [ptr] "r" (ptr),       \
			    [c0] "b" (2*192+0),     \
			   [c64] "b" (2*192+64),    \
			  [c128] "b" (2*192+128)    \
		);                            \
	} while (0)
#else
#define bgq_su3_weylnextnext_prefetch_double(addr)      \
	bgq_prefetch((char*)(addr) + 2*192 +   0);    \
	bgq_prefetch((char*)(addr) + 2*192 +  64);    \
	bgq_prefetch((char*)(addr) + 2*192 + 128)
#endif

#if BGQ_QPX
#define bgq_su3_weylnextnext_prefetch_float(addr)     \
	do {                             \
		void *ptr = (addr);          \
		asm (                        \
			"dcbt   %[c0],%[ptr]  \n" \
			"dcbt  %[c64],%[ptr]  \n" \
			"dcbt %[c128],%[ptr]  \n" \
			: : [ptr] "r" (ptr),       \
			    [c0] "b" (2*96+0),     \
			   [c64] "b" (2*96+64),    \
			  [c128] "b" (2*96+128)    \
		);                            \
	} while (0)
#else
#define bgq_su3_weylnextnext_prefetch_float(addr)      \
	bgq_prefetch((char*)(addr) + 2*96 +   0);    \
	bgq_prefetch((char*)(addr) + 2*96 +  64);    \
	bgq_prefetch((char*)(addr) + 2*96 + 128)
#endif

#if BGQ_QPX
#define bgq_su3_weyl_prefetch_float(addr)      \
	do {                             \
		void *ptr = (addr);          \
		asm (                        \
			"dcbt       0,%[ptr]  \n" \
			"dcbt  %[c64],%[ptr]  \n" \
			: : [ptr] "r" (ptr),      \
			   [c64] "b" (64)         \
		);                            \
	} while (0)
#else
#define bgq_su3_weyl_prefetch_float(addr)      \
	bgq_prefetch((char*)(addr) +   0);    \
	bgq_prefetch((char*)(addr) +  64) /* half a cacheline */
#endif

#define bgq_su3_matrix_prefetch NAME2(bgq_su3_matrix_prefetch,PRECISION)
#if BGQ_QPX
#define bgq_su3_matrix_prefetch_double(addr)     \
do { \
	void *ptr = (addr); \
	bgq_prefetch(ptr); \
	asm ( \
		"dcbt  %[c64],%[ptr]  \n" \
		"dcbt %[c128],%[ptr]  \n" \
		"dcbt %[c192],%[ptr]  \n" \
		"dcbt %[c256],%[ptr]  \n" \
		: : [ptr] "r" (ptr), \
		  [c64] "b" (64), \
		  [c128] "b" (128), \
		  [c192] "b" (192), \
		  [c256] "b" (256) \
	); \
	} while (0) /* make it consume a semicolon */
#else
#define bgq_su3_matrix_prefetch_double(addr)      \
	bgq_prefetch((char*)(addr) +   0);    \
	bgq_prefetch((char*)(addr) +  64);    \
	bgq_prefetch((char*)(addr) + 128);    \
	bgq_prefetch((char*)(addr) + 192);    \
	bgq_prefetch((char*)(addr) + 256)
#endif

#define bgq_su3_matrixnext_prefetch NAME2(bgq_su3_matrixnext_prefetch,PRECISION)
#if BGQ_QPX
#define bgq_su3_matrixnext_prefetch_double(addr)     \
	do {                             \
		void *ptr = (addr);          \
		asm (                        \
			"dcbt   %[c0],%[ptr]  \n" \
			"dcbt  %[c64],%[ptr]  \n" \
			"dcbt %[c128],%[ptr]  \n" \
			"dcbt %[c192],%[ptr]  \n" \
			"dcbt %[c256],%[ptr]  \n" \
			: : [ptr] "r" (ptr),        \
			    [c0] "b" (288+0),     \
			   [c64] "b" (288+64),    \
			  [c128] "b" (288+128),   \
			  [c192] "b" (288+192),   \
			  [c256] "b" (288+256)    \
		);                            \
	} while (0)
#else
#define bgq_su3_matrixnext_prefetch_double(addr)     \
		bgq_prefetch((char*)(addr) + 288 +   0);    \
		bgq_prefetch((char*)(addr) + 288 +  64);    \
		bgq_prefetch((char*)(addr) + 288 + 128);    \
		bgq_prefetch((char*)(addr) + 288 + 192);    \
		bgq_prefetch((char*)(addr) + 288 + 256)
#endif


#if BGQ_QPX
#define bgq_su3_matrixnext_prefetch_float(addr)     \
	do {                             \
		void *ptr = (addr);          \
		asm (                        \
			"dcbt   %[c0],%[ptr]  \n" \
			"dcbt  %[c64],%[ptr]  \n" \
			"dcbt %[c128],%[ptr]  \n" \
			: : [ptr] "r" (ptr),        \
			    [c0] "b" (144+0),     \
			   [c64] "b" (144+64),    \
			  [c128] "b" (144+128)    \
		);                            \
	} while (0)
#else
#define bgq_su3_matrixnext_prefetch_float(addr)     \
		bgq_prefetch((char*)(addr) + 144 +   0);    \
		bgq_prefetch((char*)(addr) + 144 +  64);    \
		bgq_prefetch((char*)(addr) + 144 + 128);
#endif


#define bgq_su3_matrix_prefetch_float(addr)      \
	bgq_prefetch((char*)(addr) +   0);    \
	bgq_prefetch((char*)(addr) +  64);    \
	bgq_prefetch((char*)(addr) + 128) /* quarter of a cacheline */


#define bgq_su3_spinor_zeroload NAME2(bgq_su3_spinor_zeroload,PRECISION)
#if BGQ_QPX
#define bgq_su3_spinor_zeroload_double(addr) \
	do { \
		void *ptr = (addr); \
		asm ( \
			"dcbz       0,%[ptr]  \n" \
			"dcbz  %[c64],%[ptr]  \n" \
			"dcbz %[c128],%[ptr]  \n" \
			"dcbz %[c192],%[ptr]  \n" \
			"dcbz %[c256],%[ptr]  \n" \
			"dcbz %[c320],%[ptr]  \n" \
			: : [ptr] "r" (ptr), \
			  [c64] "b" (64), \
			  [c128] "b" (128), \
			  [c192] "b" (192), \
			  [c256] "b" (256), \
			  [c320] "b" (320) \
		); \
	} while (0)
	// 384 bytes
#else
#define bgq_su3_spinor_zeroload_double(addr)     \
	do { \
		void *ptr = (addr); \
		bgq_l1_zero(ptr); \
		bgq_dcbz(ptr, 64); \
		bgq_dcbz(ptr, 128); \
		bgq_dcbz(ptr, 192); \
		bgq_dcbz(ptr, 256); \
		bgq_dcbz(ptr, 320); \
	} while (0)
	// 384 bytes
#endif


#if BGQ_QPX
#define bgq_su3_spinor_zeroload_float(addr)
#if 0
		do { \
			void *ptr = (addr); \
			bgq_l1_zero(ptr); \
			asm ( \
				"dcbz       0,%[ptr]  \n" \
				"dcbz  %[c64],%[ptr]  \n" \
				"dcbz %[c128],%[ptr]  \n" \
				: : [ptr] "r" (ptr), \
				  [c64] "b" (64), \
				  [c128] "b" (128) \
			); \
		} while (0)
		// 192 bytes
#endif
#else
#define bgq_su3_spinor_zeroload_float(addr) \
	bgq_l1_zero((char*)(addr) + 0);          \
	bgq_prefetchforwrite((char*)(addr) +  128)
	// 192
#endif



#define bgq_su3_weyl_zeroload NAME2(bgq_su3_weyl_zeroload,PRECISION)
#if BGQ_QPX
#define bgq_su3_weyl_zeroload_double(addr) \
	do { \
		void *ptr = (addr); \
		asm ( \
		"dcbz       0,%[ptr]  \n" \
		"dcbz %[c64],%[ptr]  \n" \
		"dcbz %[c128],%[ptr]  \n" \
		: : [ptr] "r" (ptr), \
		   [c64] "b" (64), \
		  [c128] "b" (128) \
		); \
	} while (0)
	// 192 bytes
#else
#define bgq_su3_weyl_zeroload_double(addr) \
	bgq_l1_zero((char*)(addr) + 0);          \
	bgq_l1_zero((char*)(addr) + 64);     \
	bgq_l1_zero((char*)(addr) + 128);
	// 192
#endif

#define bgq_su3_weyl_zeroload_float(addr)    \
	bgq_l1_zero((char*)(addr) + 0)
	// 96



#define bgq_su3_spinor_flush NAME2(bgq_su3_spinor_flush,PRECISION)
#define bgq_su3_spinor_flush_double(addr) \
	bgq_flush((char*)(addr) +   0);    \
	bgq_flush((char*)(addr) +  64);    \
	bgq_flush((char*)(addr) + 128);    \
	bgq_flush((char*)(addr) + 192);    \
	bgq_flush((char*)(addr) + 256);    \
	bgq_flush((char*)(addr) + 320)
	// 384 bytes
#define bgq_su3_spinor_flush_float(addr) \
	bgq_flush((char*)(addr) +   0);    \
	bgq_flush((char*)(addr) +  64);    \
	bgq_flush((char*)(addr) + 128)
	// 192 bytes

#define bgq_su3_weyl_flush NAME2(bgq_su3_weyl_flush,PRECISION)
#define bgq_su3_weyl_flush_double(addr) \
	bgq_flush((char*)(addr) +   0);    \
	bgq_flush((char*)(addr) +  64);    \
	bgq_flush((char*)(addr) + 128)
	// 192 bytes

#define bgq_su3_weyl_flush_float(addr)    \
	bgq_flush((char*)(addr) +   0);    \
	bgq_flush((char*)(addr) +  64)
	// 96 bytes

#define bgq_su3_matrix_flush NAME2(bgq_su3_matrix_flush,PRECISION)
#define bgq_su3_matrix_flush_double(addr) \
	bgq_flush((char*)(addr) +   0);    \
	bgq_flush((char*)(addr) +  64);    \
	bgq_flush((char*)(addr) + 128);    \
	bgq_flush((char*)(addr) + 192);    \
	bgq_flush((char*)(addr) + 256)
	// 288 bytes
#define bgq_su3_matrix_flush_float(addr) \
	bgq_flush((char*)(addr) +   0);    \
	bgq_flush((char*)(addr) +  64);    \
	bgq_flush((char*)(addr) + 128)
	// 144 bytes




#define bgq_su3_weyl_prestore NAME2(bgq_su3_weyl_prestore,PRECISION)

#ifdef BGQ_QPX
#define bgq_su3_weyl_prestore_double(addr) \
	do { \
		void *ptr = (addr); \
		asm ( \
		"dcbtst       0,%[ptr]  \n" \
		"dcbtst  %[c64],%[ptr]  \n" \
		"dcbtst %[c128],%[ptr]  \n" \
		: : [ptr] "r" (ptr), \
		   [c64] "b" (64), \
		  [c128] "b" (128) \
		); \
	} while (0)
	// 192 byte
#else
#define bgq_su3_weyl_prestore_double(addr) \
	do { \
		bgq_prefetchforwrite((char*)(addr) +   0); \
		bgq_prefetchforwrite((char*)(addr) +  64); \
		bgq_prefetchforwrite((char*)(addr) + 128); \
	while (0)
	// 192 bytes
#endif

#ifdef BGQ_QPX
#define bgq_su3_weyl_prestore_float(addr) \
	do { \
		void *ptr = (addr); \
		asm ( \
		"dcbtst       0,%[ptr]  \n" \
		"dcbtst  %[c64],%[ptr]  \n" \
		: : [ptr] "r" (ptr), \
		   [c64] "b" (64) \
		); \
	} while (0)
	// 96 byte
#else
#define bgq_su3_weyl_prestore_float(addr) \
	do { \
		bgq_prefetchforwrite((char*)(addr) +   0); \
		bgq_prefetchforwrite((char*)(addr) +  64); \
	while (0)
	// 96 bytes
#endif


#define bgq_su3_spinor_prestore NAME2(bgq_su3_spinor_prestore,PRECISION)

#ifdef BGQ_QPX
#define bgq_su3_spinor_prestore_double(addr) \
	do { \
		void *ptr = (addr); \
		asm ( \
		"dcbtst       0,%[ptr]  \n" \
		"dcbtst  %[c64],%[ptr]  \n" \
		"dcbtst %[c128],%[ptr]  \n" \
		"dcbtst %[c192],%[ptr]  \n" \
		"dcbtst %[c256],%[ptr]  \n" \
		"dcbtst %[c320],%[ptr]  \n" \
		: : [ptr] "r" (ptr), \
		   [c64] "b" (64), \
		  [c128] "b" (128), \
		  [c192] "b" (192), \
		  [c256] "b" (256), \
		  [c320] "b" (320) \
		); \
	} while (0)
	// 384 bytes
#else
#define bgq_su3_spinor_prestore_double(addr) \
	do { \
		bgq_prefetchforwrite((char*)(addr) +   0); \
		bgq_prefetchforwrite((char*)(addr) +  64); \
		bgq_prefetchforwrite((char*)(addr) + 128); \
		bgq_prefetchforwrite((char*)(addr) + 192); \
		bgq_prefetchforwrite((char*)(addr) + 256); \
		bgq_prefetchforwrite((char*)(addr) + 320); \
	while (0)
	// 384 bytes
#endif



#ifdef BGQ_QPX
#define bgq_su3_spinor_prestore_float(addr) \
	do { \
		void *ptr = (addr); \
		asm ( \
		"dcbtst       0,%[ptr]  \n" \
		"dcbtst  %[c64],%[ptr]  \n" \
		"dcbtst %[c128],%[ptr]  \n" \
		: : [ptr] "r" (ptr), \
		   [c64] "b" (64), \
		  [c128] "b" (128) \
		); \
	} while (0)
	// 384 bytes
#else
#define bgq_su3_spinor_prestore_float(addr) \
	do { \
		bgq_prefetchforwrite((char*)(addr) +   0); \
		bgq_prefetchforwrite((char*)(addr) +  64); \
		bgq_prefetchforwrite((char*)(addr) + 128); \
	while (0)
	// 384 bytes
#endif


#ifndef XLC
#include <omp.h>

static inline int Kernel_ProcessorID() {
	return omp_get_thread_num();
}

static inline int Kernel_ProcessorThreadID() {
	return 0;
}

static inline int Kernel_ProcessorCoreID() {
	return omp_get_thread_num();
}

// Memory barrier
static inline void mbar() {
	__sync_synchronize();
}

#endif


#define REORDER_BARRIER \
	__asm volatile ("");



void bgq_qpx_unittest(void);


#if BGQ_QPX
EXTERN_INLINE uint64_t bgq_wcycles() {
	return GetTimeBase();
}
#else
EXTERN_INLINE uint64_t bgq_wcycles() {
	return MPI_Wtick();
}
#endif


#define BGQ_VECTORSIZE NAME2(BGQ_VECTORSIZE,PRECISION)
#define BGQ_VECTORSIZE_double 32
#define BGQ_VECTORSIZE_float 16



#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_QPX_H_ */
