/*
 * bgq_stdreductions.c
 *
 *  Created on: Nov 10, 2012
 *      Author: meinersbur
 */

#include "bgq_qpx.h"
#include "bgq_field.h"


static inline void bgq_reduce_initkahan(bgq_params(*ks), bgq_params(*kc)) {
	bgq_zero(*ks);
	bgq_zero(*kc);
}


#define bgq_kahan_add(ks, kc, val) bgq_kahan_add_raw(bgq_vars(ks), bgq_vars(kc), bgq_vars(val))
static inline void bgq_kahan_add_raw(bgq_params(*accumulator_ks), bgq_params(*accumulator_kc), bgq_params(val)) {
	bgq_vector4double_decl(ks);
	bgq_vector4double_decl(kc);
	bgq_mov(ks, *accumulator_ks);
	bgq_mov(kc, *accumulator_kc);

	bgq_vector4double_decl(tr); // val + compensation
	bgq_vector4double_decl(ts); // next_sum
	bgq_vector4double_decl(tt); // val + next_compensation

	bgq_add(tr, val, kc); // val + carried_error
	bgq_add(ts, tr, ks); // val + carried_error + sum - error
	bgq_sub(tt, ts, ks); // val + carried_error - error
	bgq_mov(ks, ts);
	bgq_sub(kc, tr, tt); // error

    bgq_mov(*accumulator_ks, ks);
    bgq_mov(*accumulator_kc, kc);
}


static inline void bgq_reduce_prod(bgq_params(*accumulator_ks), bgq_params(*accumulator_kc), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), ucoord ic) {
	bgq_vector4double_decl(result);

	// re = a.re * b.re + a.im * b.im + re
    // im = a.re * b.im - a.im * b.re + im
	bgq_xmul(result, spinor1_v0_c0, spinor2_v0_c0);
    bgq_rxxcpnmadd(result, spinor1_v0_c0, spinor2_v0_c0, result);

    bgq_xmadd(result, spinor1_v0_c1, spinor2_v0_c1, result);
    bgq_rxxcpnmadd(result, spinor1_v0_c1, spinor2_v0_c1, result);

    bgq_xmadd(result, spinor1_v0_c2, spinor2_v0_c2, result);
    bgq_rxxcpnmadd(result, spinor1_v0_c2, spinor2_v0_c2, result);

    bgq_xmadd(result, spinor1_v1_c0, spinor2_v1_c0, result);
    bgq_rxxcpnmadd(result, spinor1_v1_c0, spinor2_v1_c0, result);

    bgq_xmadd(result, spinor1_v1_c1, spinor2_v1_c1, result);
    bgq_rxxcpnmadd(result, spinor1_v1_c1, spinor2_v1_c1, result);

    bgq_xmadd(result, spinor1_v1_c2, spinor2_v1_c2, result);
    bgq_rxxcpnmadd(result, spinor1_v1_c2, spinor2_v1_c2, result);

    bgq_xmadd(result, spinor1_v2_c0, spinor2_v2_c0, result);
    bgq_rxxcpnmadd(result, spinor1_v2_c0, spinor2_v2_c0, result);

    bgq_xmadd(result, spinor1_v2_c1, spinor2_v2_c1, result);
    bgq_rxxcpnmadd(result, spinor1_v2_c1, spinor2_v2_c1, result);

    bgq_xmadd(result, spinor1_v2_c2, spinor2_v2_c2, result);
    bgq_rxxcpnmadd(result, spinor1_v2_c2, spinor2_v2_c2, result);

    bgq_xmadd(result, spinor1_v3_c0, spinor2_v3_c0, result);
    bgq_rxxcpnmadd(result, spinor1_v3_c0, spinor2_v3_c0, result);

    bgq_xmadd(result, spinor1_v3_c1, spinor2_v3_c1, result);
    bgq_rxxcpnmadd(result, spinor1_v3_c1, spinor2_v3_c1, result);

    bgq_xmadd(result, spinor1_v3_c2, spinor2_v3_c2, result);
    bgq_rxxcpnmadd(result, spinor1_v3_c2, spinor2_v3_c2, result);

    bgq_kahan_add(accumulator_ks, accumulator_kc, result);
}


static inline void bgq_reduce_prod_r(bgq_params(*accumulator_ks), bgq_params(*accumulator_kc), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), ucoord ic) {
	bgq_vector4double_decl(result);

	bgq_mul(result, spinor1_v0_c0, spinor2_v0_c0);
	bgq_madd(result, spinor1_v0_c1, spinor2_v0_c1, result);
	bgq_madd(result, spinor1_v0_c2, spinor2_v0_c2, result);
	bgq_madd(result, spinor1_v1_c0, spinor2_v1_c0, result);
	bgq_madd(result, spinor1_v1_c1, spinor2_v1_c1, result);
	bgq_madd(result, spinor1_v1_c2, spinor2_v1_c2, result);
	bgq_madd(result, spinor1_v2_c0, spinor2_v2_c0, result);
	bgq_madd(result, spinor1_v2_c1, spinor2_v2_c1, result);
	bgq_madd(result, spinor1_v2_c2, spinor2_v2_c2, result);
	bgq_madd(result, spinor1_v3_c0, spinor2_v3_c0, result);
	bgq_madd(result, spinor1_v3_c1, spinor2_v3_c1, result);
	bgq_madd(result, spinor1_v3_c2, spinor2_v3_c2, result);

	bgq_kahan_add(accumulator_ks, accumulator_kc, result);
}


void bgq_reduce_norm(bgq_params(*ks), bgq_params(*kc), bgq_su3_spinor_params(spinor), ucoord ic) {
	bgq_vector4double_decl(siteresult);
	bgq_mul(siteresult, spinor_v0_c0, spinor_v0_c0);
	bgq_madd(siteresult, spinor_v0_c1, spinor_v0_c1, siteresult);
	bgq_madd(siteresult, spinor_v0_c2, spinor_v0_c2, siteresult);
	bgq_madd(siteresult, spinor_v1_c0, spinor_v1_c0, siteresult);
	bgq_madd(siteresult, spinor_v1_c1, spinor_v1_c1, siteresult);
	bgq_madd(siteresult, spinor_v1_c2, spinor_v1_c2, siteresult);
	bgq_madd(siteresult, spinor_v2_c0, spinor_v2_c0, siteresult);
	bgq_madd(siteresult, spinor_v2_c1, spinor_v2_c1, siteresult);
	bgq_madd(siteresult, spinor_v2_c2, spinor_v2_c2, siteresult);
	bgq_madd(siteresult, spinor_v3_c0, spinor_v3_c0, siteresult);
	bgq_madd(siteresult, spinor_v3_c1, spinor_v3_c1, siteresult);
	bgq_madd(siteresult, spinor_v3_c2, spinor_v3_c2, siteresult);

	bgq_kahan_add(ks, kc, siteresult);
}





static inline void bgq_combine_kahan(bgq_params(*ks1), bgq_params(*kc1), bgq_params(ks2), bgq_params(kc2)) {
	// We might also just sum up ks and kc together
	bgq_kahan_add(ks1, kc1, kc2);
	bgq_kahan_add(ks1, kc1, ks2);
}


#define REDUCTION_NAME bgq_spinorfield_innerprod_raw
#define REDUCTION_ARGFIELDS 2
#define REDUCTION_REDTYPES bgq_types,bgq_types
#define REDUCTION_REDARGS bgq_vars(ks),bgq_vars(kc)
#define REDUCTION_VARINIT bgq_reduce_initkahan
#define REDUCTION_SITEREDUCEFUNC bgq_reduce_prod
#define REDUCTION_COMBINEFUNC bgq_combine_kahan
#include "bgq_reduction.inc.c"

complex_double bgq_spinorfield_innerprod_local(bool isOdd, bgq_weylfield_controlblock *field1, bgq_weylfield_controlblock *field2) {
	bgq_vector4double_decl(ks);
	bgq_vector4double_decl(kc);
	bgq_spinorfield_innerprod_raw(bgq_vars(&ks), bgq_vars(&kc), isOdd, field1, field2);
	complex_double localresult = bgq_cmplxval1(kc) + bgq_cmplxval2(kc) + bgq_cmplxval1(ks) + bgq_cmplxval2(ks);
	return localresult;
}

complex_double bgq_spinorfield_innerprod_global(bool isOdd, bgq_weylfield_controlblock *field1, bgq_weylfield_controlblock *field2) {
	complex_double localresult = bgq_spinorfield_innerprod_local(isOdd, field1, field2);
	complex_double globalresult;
	MPI_Allreduce(&localresult, &globalresult, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
	return globalresult;
}


#define REDUCTION_NAME bgq_innerprod_r_raw
#define REDUCTION_ARGFIELDS 2
#define REDUCTION_REDTYPES bgq_types,bgq_types
#define REDUCTION_REDARGS bgq_vars(ks),bgq_vars(kc)
#define REDUCTION_VARINIT bgq_reduce_initkahan
#define REDUCTION_SITEREDUCEFUNC bgq_reduce_prod_r
#define REDUCTION_COMBINEFUNC bgq_combine_kahan
#include "bgq_reduction.inc.c"

double bgq_spinorfield_innerprod_r_local(bool isOdd, bgq_weylfield_controlblock *field1, bgq_weylfield_controlblock *field2) {
	bgq_vector4double_decl(ks);
	bgq_vector4double_decl(kc);
	bgq_innerprod_r_raw(bgq_vars(&ks), bgq_vars(&kc), isOdd, field1, field2);
	double localresult = bgq_elem0(kc) + bgq_elem1(kc) + bgq_elem2(kc) + bgq_elem3(kc); // This is probably useless
	localresult += bgq_elem0(ks);
	localresult += bgq_elem1(ks);
	localresult += bgq_elem2(ks);
	localresult += bgq_elem3(ks);
	return localresult;
}

double bgq_spinorfield_innerprod_r_global(bool isOdd, bgq_weylfield_controlblock *field1, bgq_weylfield_controlblock *field2) {
	double localresult = bgq_spinorfield_innerprod_r_local(isOdd, field1, field2);
	double globalresult;
	MPI_Allreduce(&localresult, &globalresult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return globalresult;
}




#define REDUCTION_NAME bgq_spinorfield_sqrnorm_raw
#define REDUCTION_ARGFIELDS 1
#define REDUCTION_REDTYPES bgq_types,bgq_types
#define REDUCTION_REDARGS bgq_vars(ks),bgq_vars(kc)
#define REDUCTION_VARINIT bgq_reduce_initkahan
#define REDUCTION_SITEREDUCEFUNC bgq_reduce_norm
#define REDUCTION_COMBINEFUNC bgq_combine_kahan
#include "bgq_reduction.inc.c"

double bgq_spinorfield_sqrnorm_local(tristate isOdd, bgq_weylfield_controlblock *field) {
	bgq_vector4double_decl(ks);
	bgq_vector4double_decl(kc);
	bgq_spinorfield_sqrnorm_raw(bgq_vars(&ks), bgq_vars(&kc), isOdd, field);
	double localresult = bgq_elem0(kc) + bgq_elem1(kc) + bgq_elem2(kc) + bgq_elem3(kc); // This is probably useless, because too small compared to ks
	localresult += bgq_elem0(ks);
	localresult += bgq_elem1(ks);
	localresult += bgq_elem2(ks);
	localresult += bgq_elem3(ks);
	return localresult;
}

double bgq_spinorfield_sqrnorm_global(tristate isOdd, bgq_weylfield_controlblock *field) {
	double localresult = bgq_spinorfield_sqrnorm_local(isOdd, field);
	double globalresult;
	MPI_Allreduce(&localresult, &globalresult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return globalresult;
}





