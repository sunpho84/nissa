/*
 * bgq_utils.h
 *
 *  Created on: Aug 4, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_UTILS_H_
#define BGQ_UTILS_H_

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#if __GNUC__ && !__GNUC_STDC_INLINE__
#define EXTERN_INLINE_DECLARATION extern inline
#define EXTERN_INLINE_DEFINITION inline
#else
#define EXTERN_INLINE_DECLARATION inline
#define EXTERN_INLINE_DEFINITION extern inline
#endif

#ifndef BGQ_UTILS_C_
#define EXTERN_INLINE EXTERN_INLINE_DECLARATION
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE EXTERN_INLINE_DEFINITION
#define EXTERN_FIELD
#define EXTERN_INIT(val) = (val)
#endif


#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic error "-Wimplicit-int"
#pragma GCC diagnostic error "-Wimplicit-function-declaration"
#endif


typedef _Complex double complexdouble;
typedef _Complex float complexfloat;

typedef _Complex double complex_double;
typedef _Complex float complex_float;

#define _CONCAT(X,Y) X##Y
#define CONCAT(X,Y) _CONCAT(X,Y)

#define CONCAT2(s1,s2) CONCAT(s1,s2)
#define CONCAT3(s1,s2,s3) CONCAT(CONCAT2(s1,s2),s3)
#define CONCAT4(s1,s2,s3,s4) CONCAT(CONCAT3(s1,s2,s3),s4)
#define CONCAT5(s1,s2,s3,s4,s5) CONCAT(s1,CONCAT4(s2,s3,s4,s5))

#define NAME2(s1,s2) CONCAT3(s1,_,s2)
#define NAME3(s1,s2,s3) NAME2(CONCAT3(s1,_,s2),s3)
#define NAME4(s1,s2,s3,s4) NAME3(CONCAT3(s1,_,s2),s3,s4)
#define NAME5(s1,s2,s3,s4,s5) NAME4(CONCAT3(s1,_,s2),s3,s4,s5)

#define MAKENAME CONCAT2(tmp,__LINE__) /* TODO: Use __COUNTER__ if available */
#define MAKENAME1(s1) CONCAT3(MAKENAME, _, s1)
#define MAKENAME2(s1,s2) CONCAT3(MAKENAME1(s1), _, s2)
#define MAKENAME3(s1,s2,s3) CONCAT3(MAKENAME2(s1,s2), _, s3)
#define MAKENAME4(s1,s2,s3,s4) CONCAT3(MAKENAME3(s1,s2,s3), _, s4)
#define MAKENAME5(s1,s2,s3,s4,s5) CONCAT3(MAKENAME4(s1,s2,s3,s4), _, s5)
#define MAKENAME6(s1,s2,s3,s4,s5,s6) CONCAT3(MAKENAME5(s1,s2,s3,s4,s5), _, s6)

#ifndef STRINGIFY
#define STRINGIFY(V) #V
#endif
#ifndef TOSTRING
#define TOSTRING(V) STRINGIFY(V)
#endif




#define MPI_CHECK(RTNCODE)                                                                 \
	do {                                                                                   \
    	int mpi_rtncode = (RTNCODE);                                                       \
    	if (mpi_rtncode != MPI_SUCCESS) {                                                  \
			fprintf(stderr, "MPI call %s at %s:%d failed: errorcode %d\n", #RTNCODE, __FILE__, __LINE__, mpi_rtncode);  \
			assert(!"MPI call " #RTNCODE " failed");                                       \
			abort();                                                                       \
		}                                                                                  \
	} while (0)

#define L1P_CHECK(RTNCODE)                                                                 \
	do {                                                                                   \
    	int mpi_rtncode = (RTNCODE);                                                       \
    	if (mpi_rtncode != MPI_SUCCESS) {                                                  \
			fprintf(stderr, "L1P call %s at %s:%d failed: errorcode %d\n", #RTNCODE, __FILE__, __LINE__, mpi_rtncode);  \
			assert(!"L1P call " #RTNCODE " failed");                                       \
			abort();                                                                       \
		}                                                                                  \
	} while (0)

EXTERN_INLINE int get_MPI_count(MPI_Status *status) {
	int count;
	MPI_CHECK(MPI_Get_count(status, MPI_BYTE, &count));
	return count;
}

#define WORKLOAD_DECL(COUNTER, TOTAL) \
	size_t xyz_counter = (COUNTER);       \
	size_t xyz_isntance = xyz_counter;    \
	size_t xyz_orig;                      \
	size_t xyz_torig;                     \
	size_t xyz_total = (TOTAL);           \
	size_t xyz_param;                     \
	assert(xyz_counter >= 0);          \
	assert(xyz_counter < xyz_total)

#if 0
if ((xyz_counter < 0) || (xyz_counter > xyz_total)) {   \
	node_print("Invalid counter %d\n", xyz_counter);  \
}  \

if (xyz_counter==0) { \
	fprintf(stderr, "xyz_counter=%d xyz_total=%d xyz_isntance=%d BODY_ZLINES=%d\n", xyz_counter, xyz_total, xyz_isntance, BODY_ZLINES); \
} \
//
#endif

// true:  xyz = 0            .. TRUE_COUNT -> 0 .. TRUE_COUNT
// false: xyz = TRUE_COUNT+1 .. xyz_total  -> 0 .. xyz_total-TRUE_COUNT
//TODO: Can also do a variant without conditional
#define WORKLOAD_SPLIT(TRUE_COUNT) \
	((xyz_counter < TRUE_COUNT) ? (                    \
		xyz_total = (TRUE_COUNT),              \
		true                                   \
	) : (                                      \
		xyz_total = xyz_total - (TRUE_COUNT),  \
		xyz_counter = xyz_counter - (TRUE_COUNT),              \
		assert(xyz_total >= 0 && "Expected more to split"),   \
		false                                  \
	))
/*
 xyz_orig		xyz				return
 0				0				1
 1				1				1
 2				2				1
 ...
 TRUE_COUNT-1	TRUE_COUNT-1	1
 TRUE_COUNT		0				0
 1				0
 2				0
 ...
 xyz_torig-1						0
 */

#define WORKLOAD_PARAM(LENGTH)            \
	(assert((mod(xyz_total, (LENGTH))==0) && "Loop bounds must be a multiple of this parameter"), \
	 xyz_param = (LENGTH),                \
	 xyz_orig = xyz_counter,                      \
	 xyz_torig = xyz_total,               \
	 xyz_total = xyz_torig / xyz_param,   \
	 xyz_counter = xyz_orig / xyz_param,          \
	 (xyz_orig % xyz_param))
/*
 xyz_orig	xyz			return (=PARAM)
 0			0			0
 1			0			1
 2			0			2
 ...
 LENGTH-1	0			LENGTH-1
 LENGTH		1			0
 1			1
 1			2
 ...
 xyz_torig	xyz_total	LENGTH

 xyz_orig =	xyz*LENGTH + return
 */

#define WORKLOAD_SECTION(SECTIONS)        \
	(xyz_param = (SECTIONS),              \
	 xyz_orig  = xyz_counter,                     \
	 xyz_torig = xyz_total,               \
	 xyz_total = xyz_param,  			  \
	 xyz_counter = xyz_orig / (xyz_torig/xyz_param),          \
	 (xyz_orig % (xyz_torig/xyz_param)))
/*
 xyz_orig	xyz			return
 0			0			0
 1			0			1
 2			0			2
 ...
 0
 1			0
 1			1
 1			2
 ...
 xyz_torig-1	SECTIONS-1
 xyz_torig	SECTIONS	xyz_torig/SECTIONS
 (=xyz_total)

 xyz_orig =	xyz*(xyz_torig/SECTIONS) + return
 */

#define WORKLOAD_TILE(TILES)   					   \
		(xyz_param = (TILES),                      \
		 xyz_orig = xyz_counter,                           \
		 xyz_torig = xyz_total,                    \
		 xyz_total = TILES,     				   \
		 xyz_counter = (xyz_orig % (xyz_torig / xyz_param)), \
		 xyz_orig / (xyz_torig / xyz_param))
/*
 xyz_orig	xyz	(=TILE)	return
 0			0			0
 1			1			0
 2			2			0
 ...
 TILES-1	TILES-1		0
 TILES		0			1
 1			1
 2			1
 ...
 xyz_torig	TILES		xyz_torig/TILES
 (=xyz_total)

 xyz_orig =	xyz +		return*TILES
 */

#define WORKLOAD_CHUNK(LENGTH)                 \
	(assert((mod(xyz_total, (LENGTH))==0) && "Loop bounds must be a multiple of this parameter"), \
	 xyz_param = (LENGTH),                     \
	 xyz_orig = xyz_counter,                           \
	 xyz_torig = xyz_total,                    \
	 xyz_total = xyz_torig / xyz_param,        \
	 xyz_counter = (xyz_orig % xyz_total),               \
	 xyz_orig / xyz_total)
/*
 xyz_orig	xyz			return (=CHUNK)
 0			0			0
 1			1			0
 2			2			0
 ...
 0
 0			1
 1			1
 2			1
 ...
 xyz_torig-1	xyz_total-1	LENGTH-1
 xyz_torig	xyz_total	LENGTH
 (=xyz_torig/LENGTH)

 xyz_orig =	xyz + 		return*xyz_total
 */

#ifdef NDEBUG
#define WORKLOAD_CHECK
#else
#define WORKLOAD_CHECK \
	if (xyz_counter!=0 && g_proc_id==0) { \
		fprintf(stderr, "xyz_counter=%zu xyz_total=%zu xyz_isntance=%zu\n", xyz_counter, xyz_total, xyz_isntance); \
	} \
	assert(xyz_counter==0); \
	assert(xyz_total==1);
#endif

#ifndef BGQMOD
#define BGQMOD 0
#endif

// dividend >= -1
// divisor > 0
// 0 <= result < divisor
EXTERN_INLINE int moddown(const int dividend, const int divisor) {
#if BGQMOD==0
	// Compilers can therefore optimize it to bit-operations specific to the target machine
	// This is not possible with the %-operator on signed operands because it required the result to have the sign of the dividend (hence its the remainder, not a modulus)
	unsigned int udividend = dividend + divisor; // In all our use cases, dividend is at lowest -1, so just add one divisor such that result stays the same but everything is positive
	unsigned int udivisor = divisor;
	return udividend % udivisor; // Use unsigned operation here to enable the optimizer to use bit-tricks (like dividend&1 if divisor==2)
#elif BGQMOD==1
	int result = dividend % divisor;
	if (result < 0)
		result += divisor;
	return result;
#elif BGQMOD==2
	if (dividend >= 0)
		return dividend % divisor;
	else
		return dividend % divisor + divisor;
#elif BGQMOD==3
	return ((dividend % divisor) + divisor) % divisor;
#endif
}

 // dividend >= 0
 // divisor > 0
 // 0 <= result < divisor
EXTERN_INLINE int mod(const int dividend, const int divisor) {
  	unsigned int udividend = dividend;
 	unsigned int udivisor = divisor;
 	return udividend % udivisor; // Use unsigned operation here to enable the optimizer to use bit-tricks (like dividend&1 if divisor==2)
}


#ifndef BGQDIV
#define BGQDIV 0
#endif

// Always round towards negative infinity
// dividend >= -1
// divisor > 0
// result >= -1
EXTERN_INLINE int divdown(const int dividend, const int divisor) {
#if BGQDIV==0
	unsigned int udividend = dividend + divisor; // Dividend can be -1
	unsigned int udivisor = divisor;
	return (int)(udividend / udivisor) - 1;
#elif BGQDIV==1
	return (dividend - mod(dividend,divisor)) / divisor;
#elif BGQDIV==2
	if (dividend >= 0)
		return dividend / divisor;
	else
		return (dividend - divisor + 1) / divisor;
#endif
}



#define BGQ_ENTER_FUNC                                     \
	if (g_proc_id==0) {                                    \
		fprintf(stderr, "MK ENTER_FUNC %s\n",  __func__);  \
	}

#define BGQ_BEACON \
	{  \
	static bool CONCAT(beacon, __LINE__); \
	if (!CONCAT(beacon, __LINE__)) \
		master_print("BEACON: File: %s, line %d, func: %s\n", __FILE__, __LINE__, __func__); \
	CONCAT(beacon, __LINE__) = true; \
	}


#define master_print(...)              \
	if (g_proc_id == 0)                 \
		if (omp_get_thread_num() == 0)  \
			fprintf(stderr, __VA_ARGS__)

#define node_print(...)                          \
	do {                                          \
		fprintf(stderr, "rank=%d thread=%d ", g_proc_id, omp_get_thread_num()); \
		fprintf(stderr, __VA_ARGS__);  \
	} while (0)

#define master_error(errcode, ...) \
	do {                            \
		master_print(__VA_ARGS__);  \
		assert(false);              \
		exit(errcode);              \
	} while (0)

#define lengthof(X) (sizeof(X)/sizeof((X)[0]))

void *malloc_aligned(size_t size, size_t alignment);


EXTERN_INLINE double max_double(double lhs, double rhs) {
	if (lhs > rhs)
		return lhs;
	return rhs;
}


EXTERN_INLINE double min_double(double lhs, double rhs) {
	if (lhs > rhs)
		return rhs;
	return lhs;
}


EXTERN_INLINE size_t min_sizet(size_t lhs, size_t rhs) {
	if (lhs > rhs)
		return rhs;
	return lhs;
}


EXTERN_INLINE size_t max_sizet(size_t lhs, size_t rhs) {
	if (lhs > rhs)
		return lhs;
	return rhs;
}


// from stackoverflow
// log_2(_v)
EXTERN_INLINE int ilog(unsigned int _v) {
	int ret;
	int m;
	ret=!!_v;
	m=!!(_v&0xFFFF0000)<<4;
	_v>>=m;
	ret|=m;
	m=!!(_v&0xFF00)<<3;
	_v>>=m;
	ret|=m;
	m=!!(_v&0xF0)<<2;
	_v>>=m;
	ret|=m;
	m=!!(_v&0xC)<<1;
	_v>>=m;
	ret|=m;
	ret+=!!(_v&0x2);
	return ret;
}


void opaque_func_call();

EXTERN_INLINE double sqr(double val) {
	return val*val;
}


#define ADD_PADDING(addr,alignment) (void*)(((uintptr_t)(addr) + ((uintptr_t)(alignment)-1)) & ~((uintptr_t)(alignment)-1))


EXTERN_INLINE size_t gcd_sizet(size_t a, size_t b) {
	assert(1 <= a);
	assert(1 <= b);
	while (b) {
		size_t t = b;
		b = a % b;
		a = t;
	}
	return a;
}


EXTERN_INLINE size_t lcm_sizet(size_t a, size_t b) {
	assert(1 <= a);
	assert(1 <= b);
	assert((a*b) % gcd_sizet(a,b) == 0);
	return (a*b) /  gcd_sizet(a,b);
}

#ifdef __GNUC__
#define UNREACHABLE \
	{ \
	assert(!"Unreachable"); \
	__builtin_unreachable(); \
	abort(); \
	}
#else
#define UNREACHABLE \
	{ \
	assert(!"Unreachable"); \
	abort(); \
	}
#endif


#define REMOVE_PAREN(...) __VA_ARGS__
#define REMOVE_PAREN2(...) REMOVE_PAREN __VA_ARGS__

#define _NUM_ARGS2(X,X64,X63,X62,X61,X60,X59,X58,X57,X56,X55,X54,X53,X52,X51,X50,X49,X48,X47,X46,X45,X44,X43,X42,X41,X40,X39,X38,X37,X36,X35,X34,X33,X32,X31,X30,X29,X28,X27,X26,X25,X24,X23,X22,X21,X20,X19,X18,X17,X16,X15,X14,X13,X12,X11,X10,X9,X8,X7,X6,X5,X4,X3,X2,X1,N,...) N
#define NUM_ARGS(...) _NUM_ARGS2(0, __VA_ARGS__ ,64,63,62,61,60,59,58,57,56,55,54,53,52,51,50,49,48,47,46,45,44,43,42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0)


// Thanks to Jens Gustedt
// Taken from http://gustedt.wordpress.com/2010/06/08/detect-empty-macro-arguments/
#define _ARG16(_0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, ...) _15
#define HAS_COMMA(...) _ARG16(__VA_ARGS__, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, more than 16 arguments used/*dummy to conform to C99*/)
#define _TRIGGER_PARENTHESIS_(...) ,
#define ISEMPTY(...)                                                    \
_ISEMPTY(                                                               \
          /* test if there is just one argument, eventually an empty    \
             one */                                                     \
          HAS_COMMA(__VA_ARGS__),                                       \
          /* test if _TRIGGER_PARENTHESIS_ together with the argument   \
             adds a comma */                                            \
          HAS_COMMA(_TRIGGER_PARENTHESIS_ __VA_ARGS__),                 \
          /* test if the argument together with a parenthesis           \
             adds a comma */                                            \
          HAS_COMMA(__VA_ARGS__ (/*empty*/)),                           \
          /* test if placing it between _TRIGGER_PARENTHESIS_ and the   \
             parenthesis adds a comma */                                \
          HAS_COMMA(_TRIGGER_PARENTHESIS_ __VA_ARGS__ (/*empty*/))      \
          )
#define PASTE5(_0, _1, _2, _3, _4) _0 ## _1 ## _2 ## _3 ## _4
#define _ISEMPTY(_0, _1, _2, _3) HAS_COMMA(PASTE5(_IS_EMPTY_CASE_, _0, _1, _2, _3))
#define _IS_EMPTY_CASE_0001 ,


#define _UNPARENTHESIZE0(...) __VA_ARGS__
#define _UNPARENTHESIZE1(...) REMOVE_PAREN __VA_ARGS__
#define UNPARENTHESIZE(...) CONCAT(_UNPARENTHESIZE, HAS_COMMA(_TRIGGER_PARENTHESIS_ __VA_ARGS__))(__VA_ARGS__)


#define _IFELSE_0(ifcase,elsecase) elsecase
#define _IFELSE_1(ifcase,elsecase) ifcase
#define IF_EMPTY(ifcase,elsecase,...) CONCAT(_IFELSE_, ISEMPTY(__VA_ARGS__)) (ifcase,elsecase)


#define _IFELSE_INPAREN_0(ifcase,elsecase) REMOVE_PAREN elsecase
#define _IFELSE_INPAREN_1(ifcase,elsecase) REMOVE_PAREN ifcase
#define IF_EMPTY_INPAREN(ifcase,elsecase, ...) CONCAT(_IFELSE_INPAREN_, ISEMPTY(__VA_ARGS__)) (ifcase,elsecase)
#define IF_TWOARGSORMORE(ifcase,elsecase,...) CONCAT(_IFELSE_, HAS_COMMA(__VA_ARGS__)) (ifcase,elsecase)


// True recursion is impossible with preprocessor
//TODO: Version that accepts comma in parenthesis in arguments
#define _PREPROCESSOR_FOREACH_VARARGS0(...) /* dummy version, never used */
#define _PREPROCESSOR_FOREACH_VARARGS1(infix,callback,one) callback(one)
#define _PREPROCESSOR_FOREACH_VARARGS2(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS1(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS3(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS2(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS4(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS3(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS5(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS4(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS6(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS5(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS7(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS6(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS8(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS7(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS9(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS8(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS10(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS9(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS11(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS10(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS12(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS11(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS13(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS12(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS14(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS13(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS15(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS14(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS16(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS15(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS17(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS16(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS18(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS17(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS19(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS18(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS20(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS19(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS21(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS20(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS22(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS21(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS23(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS22(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS24(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS23(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS25(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS24(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS26(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS25(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS27(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS26(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS28(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS27(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS29(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS28(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS30(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS29(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS31(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS30(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS32(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS31(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS33(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS32(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS34(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS33(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS35(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS34(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS36(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS35(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS37(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS36(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS38(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS37(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS39(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS38(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS40(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS39(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS41(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS40(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS42(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS41(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS43(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS42(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS44(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS43(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS45(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS44(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS46(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS45(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS47(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS46(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS48(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS47(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS49(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS48(infix,callback,__VA_ARGS__)
#define _PREPROCESSOR_FOREACH_VARARGS50(infix,callback,one,...) callback(one) UNPARENTHESIZE(infix) _PREPROCESSOR_FOREACH_VARARGS49(infix,callback,__VA_ARGS__)

#define PREPROCESSOR_FOREACH(prefix,infix,postfix,emptycase,callback, ...) REMOVE_PAREN2(IF_EMPTY((UNPARENTHESIZE(emptycase)),(UNPARENTHESIZE(prefix) CONCAT(_PREPROCESSOR_FOREACH_VARARGS, NUM_ARGS(__VA_ARGS__))(infix,callback,__VA_ARGS__) UNPARENTHESIZE(postfix)),__VA_ARGS__))

#define IDENTITY(...) __VA_ARGS__
#define EMPTY(...)

#define PP_UNARG0(one, ...) one
#define PP_UNARG1(x1,one, ...) one
#define PP_UNARG2(x1,x2,one, ...) one
#define PP_UNARG3(x1,x2,x3,one, ...) one
#define PP_UNARG4(x1,x2,x3,x4,one, ...) one
#define PP_UNARG5(x1,x2,x3,x4,x5,one, ...) one
#define PP_UNARG6(x1,x2,x3,x4,x5,x6,one, ...) one
#define PP_UNARG7(x1,x2,x3,x4,x5,x6,x7,one, ...) one
#define PP_UNARG8(x1,x2,x3,x4,x5,x6,x7,x8,one, ...) one
#define PP_UNARG9(x1,x2,x3,x4,x5,x6,x7,x8,x9,one, ...) one

#define PP_GETARG0(...) PP_UNARG0(__VA_ARGS__,/*dummy*/)
#define PP_GETARG1(...) PP_UNARG1(__VA_ARGS__,/*dummy*/)
#define PP_GETARG2(...) PP_UNARG2(__VA_ARGS__,/*dummy*/)
#define PP_GETARG3(...) PP_UNARG3(__VA_ARGS__,/*dummy*/)
#define PP_GETARG4(...) PP_UNARG4(__VA_ARGS__,/*dummy*/)
#define PP_GETARG5(...) PP_UNARG5(__VA_ARGS__,/*dummy*/)
#define PP_GETARG6(...) PP_UNARG6(__VA_ARGS__,/*dummy*/)
#define PP_GETARG7(...) PP_UNARG7(__VA_ARGS__,/*dummy*/)
#define PP_GETARG8(...) PP_UNARG8(__VA_ARGS__,/*dummy*/)
#define PP_GETARG9(...) PP_UNARG9(__VA_ARGS__,/*dummy*/)

#define PP_GETARG(arg, ...) CONCAT(PP_GETARG,arg) (__VA_ARGS__,/*dummy*/)

#define PP_ZIP0(macro, infix)
#define PP_ZIP1(macro, infix, list1,list2) macro(PP_GETARG0(REMOVE_PAREN list1),PP_GETARG0(REMOVE_PAREN list2))
#define PP_ZIP2(macro, infix, list1,list2) PP_ZIP1(macro,infix,list1,list2) UNPARENTHESIZE(infix) macro(PP_GETARG1(REMOVE_PAREN list1),PP_GETARG1(REMOVE_PAREN list2))
#define PP_ZIP3(macro, infix, list1,list2) PP_ZIP2(macro,infix,list1,list2) UNPARENTHESIZE(infix) macro(PP_GETARG2(REMOVE_PAREN list1),PP_GETARG2(REMOVE_PAREN list2))
#define PP_ZIP4(macro, infix, list1,list2) PP_ZIP3(macro,infix,list1,list2) UNPARENTHESIZE(infix) macro(PP_GETARG3(REMOVE_PAREN list1),PP_GETARG3(REMOVE_PAREN list2))
#define PP_ZIP5(macro, infix, list1,list2) PP_ZIP4(macro,infix,list1,list2) UNPARENTHESIZE(infix) macro(PP_GETARG4(REMOVE_PAREN list1),PP_GETARG4(REMOVE_PAREN list2))
#define PP_ZIP6(macro, infix, list1,list2) PP_ZIP5(macro,infix,list1,list2) UNPARENTHESIZE(infix) macro(PP_GETARG5(REMOVE_PAREN list1),PP_GETARG5(REMOVE_PAREN list2))
#define PP_ZIP7(macro, infix, list1,list2) PP_ZIP6(macro,infix,list1,list2) UNPARENTHESIZE(infix) macro(PP_GETARG6(REMOVE_PAREN list1),PP_GETARG6(REMOVE_PAREN list2))
#define PP_ZIP8(macro, infix, list1,list2) PP_ZIP7(macro,infix,list1,list2) UNPARENTHESIZE(infix) macro(PP_GETARG7(REMOVE_PAREN list1),PP_GETARG7(REMOVE_PAREN list2))
#define PP_ZIP9(macro, infix, list1,list2) PP_ZIP8(macro,infix,list1,list2) UNPARENTHESIZE(infix) macro(PP_GETARG8(REMOVE_PAREN list1),PP_GETARG8(REMOVE_PAREN list2))
#define PP_ZIP10(macro, infix, list1,list2) PP_ZIP9(macro,infix,list1,list2) UNPARENTHESIZE(infix) macro(PP_GETARG9(REMOVE_PAREN list1),PP_GETARG9(REMOVE_PAREN list2))
#define PP_ZIP(macro, infix, list1, list2) CONCAT(PP_ZIP,NUM_ARGS(REMOVE_PAREN list1))(macro,infix,list1,list2)





#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_UTILS_H_ */
