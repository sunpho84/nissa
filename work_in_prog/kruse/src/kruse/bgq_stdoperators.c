/*
 * bgq_stdoperators.c
 *
 *  Created on: Nov 21, 2012
 *      Author: meinersbur
 */

#include "bgq_stdoperators.h"

#include "bgq_field.h"
#include "bgq_qpx.h"






static inline void bgq_site_rmul_plain_sub(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), bgq_params(qr), ucoord ic) {
	bgq_su3_spinor_decl(result);
	bgq_su3_spinor_decl(rmul1);

	// Unfortunately, there is no fused-multiply-sub for this
	//NOTE: there is, and the compiler combines them automatically!
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


static inline void bgq_site_set(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_su3_spinor_mov(result, spinor);

	bgq_su3_spinor_mov(*target, result);
}


static inline void bgq_site_rmul(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor), bgq_params(qc), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_mul(result_v0_c0, qc, spinor_v0_c0);
	bgq_mul(result_v0_c1, qc, spinor_v0_c1);
	bgq_mul(result_v0_c2, qc, spinor_v0_c2);
	bgq_mul(result_v1_c0, qc, spinor_v1_c0);
	bgq_mul(result_v1_c1, qc, spinor_v1_c1);
	bgq_mul(result_v1_c2, qc, spinor_v1_c2);
	bgq_mul(result_v2_c0, qc, spinor_v2_c0);
	bgq_mul(result_v2_c1, qc, spinor_v2_c1);
	bgq_mul(result_v2_c2, qc, spinor_v2_c2);
	bgq_mul(result_v3_c0, qc, spinor_v3_c0);
	bgq_mul(result_v3_c1, qc, spinor_v3_c1);
	bgq_mul(result_v3_c2, qc, spinor_v3_c2);

	bgq_su3_spinor_mov(*target, result);
}


static inline void bgq_site_rmul_plain_add(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), bgq_params(qr), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_madd(result_v0_c0, qr, spinor1_v0_c0, spinor2_v0_c0);
	bgq_madd(result_v0_c1, qr, spinor1_v0_c1, spinor2_v0_c1);
	bgq_madd(result_v0_c2, qr, spinor1_v0_c2, spinor2_v0_c2);
	bgq_madd(result_v1_c0, qr, spinor1_v1_c0, spinor2_v1_c0);
	bgq_madd(result_v1_c1, qr, spinor1_v1_c1, spinor2_v1_c1);
	bgq_madd(result_v1_c2, qr, spinor1_v1_c2, spinor2_v1_c2);
	bgq_madd(result_v2_c0, qr, spinor1_v2_c0, spinor2_v2_c0);
	bgq_madd(result_v2_c1, qr, spinor1_v2_c1, spinor2_v2_c1);
	bgq_madd(result_v2_c2, qr, spinor1_v2_c2, spinor2_v2_c2);
	bgq_madd(result_v3_c0, qr, spinor1_v3_c0, spinor2_v3_c0);
	bgq_madd(result_v3_c1, qr, spinor1_v3_c1, spinor2_v3_c1);
	bgq_madd(result_v3_c2, qr, spinor1_v3_c2, spinor2_v3_c2);

	bgq_su3_spinor_mov(*target, result);
}


static inline void bgq_site_gamma5(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_mov(result_v0_c0, spinor_v0_c0);
	bgq_mov(result_v0_c1, spinor_v0_c1);
	bgq_mov(result_v0_c2, spinor_v0_c2);
	bgq_mov(result_v1_c0, spinor_v1_c0);
	bgq_mov(result_v1_c1, spinor_v1_c1);
	bgq_mov(result_v1_c2, spinor_v1_c2);
	bgq_neg(result_v2_c0, spinor_v2_c0);
	bgq_neg(result_v2_c1, spinor_v2_c1);
	bgq_neg(result_v2_c2, spinor_v2_c2);
	bgq_neg(result_v3_c0, spinor_v3_c0);
	bgq_neg(result_v3_c1, spinor_v3_c1);
	bgq_neg(result_v3_c2, spinor_v3_c2);

	bgq_su3_spinor_mov(*target, result);
}


static inline void bgq_site_add(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_su3_vadd(result_v0, spinor1_v0, spinor2_v0);
	bgq_su3_vadd(result_v1, spinor1_v1, spinor2_v1);
	bgq_su3_vadd(result_v2, spinor1_v2, spinor2_v2);
	bgq_su3_vadd(result_v3, spinor1_v3, spinor2_v3);

	bgq_su3_spinor_mov(*target, result);
}

#define OPERATOR_NAME bgq_spinorfield_add
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_add
#include "bgq_operator.inc.c"


static inline void bgq_site_sub(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_su3_vsub(result_v0, spinor1_v0, spinor2_v0);
	bgq_su3_vsub(result_v1, spinor1_v1, spinor2_v1);
	bgq_su3_vsub(result_v2, spinor1_v2, spinor2_v2);
	bgq_su3_vsub(result_v3, spinor1_v3, spinor2_v3);

	bgq_su3_spinor_mov(*target, result);
}

#define OPERATOR_NAME bgq_spinorfield_sub
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_sub
#include "bgq_operator.inc.c"



#define OPERATOR_NAME bgq_spinorfield_set
#define OPERATOR_ARGFIELDS 0
#define OPERATOR_VECSITEFUNC bgq_site_set
#define OPERATOR_EXTRAPARMS bgq_su3_spinor_params(spinor)
#define OPERATOR_EXTRAARGS bgq_su3_spinor_vars(spinor)
#include "bgq_operator.inc.c"


#define OPERATOR_NAME bgq_spinorfield_copy_raw
#define OPERATOR_ARGFIELDS 1
#define OPERATOR_VECSITEFUNC bgq_site_set
#include "bgq_operator.inc.c"

void bgq_spinorfield_copy(bgq_weylfield_controlblock *target, bgq_spinorfield_layout targetLayout, bgq_weylfield_controlblock *source) {
	assert(target);
	assert(source);
	assert(source != target);

	switch (targetLayout) {
	case ly_full_double:
		if (source->has_fulllayout_double) {
			bgq_spinorfield_prepareRead(source, tri_unknown, false, true, false, false, false);
			bgq_spinorfield_prepareWrite(target, source->isOdd, ly_full_double, false);
			bgq_master_memcpy(target->sec_fullspinor_double, source->sec_fullspinor_double, PHYSICAL_VOLUME * sizeof(*source->sec_fullspinor_double));
		} else {
			bgq_spinorfield_copy_raw_double(target, source->isOdd, source);
		}
		break;
	case ly_full_float:
		if (source->has_fulllayout_float) {
			bgq_spinorfield_prepareRead(source, tri_unknown, false, false, true, false, false);
			bgq_spinorfield_prepareWrite(target, source->isOdd, ly_full_float, false);
			bgq_master_memcpy(target->sec_fullspinor_float, source->sec_fullspinor_float, PHYSICAL_VOLUME * sizeof(*source->sec_fullspinor_float));
		} else {
			bgq_spinorfield_copy_raw_float(target, source->isOdd, source);
		}
		break;
	case ly_legacy:
		bgq_spinorfield_prepareRead(source, false, false, false, false, false, true);
		bgq_spinorfield_prepareWrite(target, source->isOdd, ly_legacy, false);
		bgq_master_memcpy(target->legacy_field, source->legacy_field, VOLUME/2 * sizeof(*source->legacy_field));
		break;
	default:
		master_error(1, "Cannot copy to layout");
	}

	assert(bgq_spinorfield_prepareRead(target, source->isOdd, true, true, true, true, true) == targetLayout);
}



static inline void bgq_site_imul(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor), bgq_params(qz), bgq_params(qw), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_su3_cvmul(result_v0,qz,spinor_v0);
	bgq_su3_cvmul(result_v1,qz,spinor_v1);
	bgq_su3_cvmul(result_v2,qw,spinor_v2);
	bgq_su3_cvmul(result_v3,qw,spinor_v3);

	bgq_su3_spinor_mov(*target, result);
}

#define OPERATOR_NAME bgq_spinorfield_imul_raw
#define OPERATOR_ARGFIELDS 1
#define OPERATOR_VECSITEFUNC bgq_site_imul
#define OPERATOR_EXTRAPARMS bgq_params(qz),bgq_params(qw)
#define OPERATOR_EXTRAARGS bgq_vars(qz),bgq_vars(qw)
#include "bgq_operator.inc.c"

void bgq_spinorfield_imul_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield, complex_double z, complex_double w) {
	bgq_vector4double_decl(qz);
	bgq_cconst(qz, creal(z), cimag(z));
	bgq_vector4double_decl(qw);
	bgq_cconst(qw, creal(w), cimag(w));
	bgq_spinorfield_imul_raw_double(targetfield, isOdd, sourcefield, bgq_vars(qz), bgq_vars(qw));
}

void bgq_spinorfield_imul_float(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield, complex_double z, complex_double w) {
	bgq_vector4double_decl(qz);
	bgq_cconst(qz, creal(z), cimag(z));
	bgq_vector4double_decl(qw);
	bgq_cconst(qw, creal(w), cimag(w));
	bgq_spinorfield_imul_raw_float(targetfield, isOdd, sourcefield, bgq_vars(qz), bgq_vars(qw));
}



#define OPERATOR_NAME bgq_spinorfield_rmul_plain_add_raw
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_rmul_plain_add
#define OPERATOR_EXTRAPARMS bgq_params(qr)
#define OPERATOR_EXTRAARGS bgq_vars(qr)
#include "bgq_operator.inc.c"


void bgq_spinorfield_rmul_plain_add_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *spinorfield1, bgq_weylfield_controlblock *spinorfield2, double r) {
	bgq_vector4double_decl(qr);
	bgq_cconst(qr, r, r);
	bgq_spinorfield_rmul_plain_add_raw_double(targetfield, isOdd, spinorfield1, spinorfield2, bgq_vars(qr));
}





#define OPERATOR_NAME bgq_spinorfield_gamma5
#define OPERATOR_ARGFIELDS 1
#define OPERATOR_VECSITEFUNC bgq_site_gamma5
#include "bgq_operator.inc.c"





#define OPERATOR_NAME bgq_spinorfield_rmul_raw
#define OPERATOR_ARGFIELDS 1
#define OPERATOR_VECSITEFUNC bgq_site_rmul
#define OPERATOR_EXTRAPARMS bgq_params(qr)
#define OPERATOR_EXTRAARGS bgq_vars(qr)
#include "bgq_operator.inc.c"

void bgq_spinorfield_rmul_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield, double r) {
	bgq_vector4double_decl(qr);
	bgq_cconst(qr, r, r);
	bgq_spinorfield_rmul_raw_double(targetfield, isOdd, sourcefield, bgq_vars(qr));
}

void bgq_spinorfield_rmul_float(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield, double r) {
	bgq_vector4double_decl(qr);
	bgq_cconst(qr, r, r);
	bgq_spinorfield_rmul_raw_float(targetfield, isOdd, sourcefield, bgq_vars(qr));
}



#define OPERATOR_NAME bgq_spinorfield_rmul_plain_sub_raw
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_rmul_plain_sub
#define OPERATOR_EXTRAPARMS bgq_params(qr)
#define OPERATOR_EXTRAARGS bgq_vars(qr)
#include "bgq_operator.inc.c"

void bgq_spinorfield_rmul_plain_sub_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield1, bgq_weylfield_controlblock *sourcefield2, double r) {
	bgq_vector4double_decl(qr);
	bgq_cconst(qr, r, r);
	bgq_spinorfield_rmul_plain_sub_raw_double(targetfield, isOdd, sourcefield1, sourcefield2, bgq_vars(qr));
}

void bgq_spinorfield_rmul_plain_sub_float(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield1, bgq_weylfield_controlblock *sourcefield2, double r) {
	bgq_vector4double_decl(qr);
	bgq_cconst(qr, r, r);
	bgq_spinorfield_rmul_plain_sub_raw_float(targetfield, isOdd, sourcefield1, sourcefield2, bgq_vars(qr));
}




static inline void bgq_site_cmul_plain_sub_gamma5(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), bgq_params(qc), bgq_params(qw), ucoord ic) {
	bgq_su3_spinor_decl(result);
	bgq_su3_spinor_decl(imul1);

	assert(bgq_cmplxval1(qc) == bgq_cmplxval2(qc));

	bgq_su3_cvmul(imul1_v0, qc, spinor1_v0);
	bgq_su3_cvmul(imul1_v1, qc, spinor1_v1);
	bgq_su3_vsub(result_v0, imul1_v0, spinor2_v0);
	bgq_su3_vsub(result_v1, imul1_v1, spinor2_v1);

#if 1
	assert(bgq_cmplxval1(qw) == bgq_cmplxval2(qw));
	assert(bgq_cmplxval1(qc) == -bgq_cmplxval1(qw));

// result.re = spinor2.re + qw.re * spinor1.re + qw.im * spinor1.im
// result.im = spinor2.im + qw.re * spinor1.im - qw.im * spinor1.re

	bgq_xmadd(result_v2_c0, qw, spinor1_v2_c0, spinor2_v2_c0);
	bgq_xxcpnmadd(result_v2_c0, spinor1_v2_c0, qw, result_v2_c0);
	bgq_xmadd(result_v2_c1, qw, spinor1_v2_c1, spinor2_v2_c1);
	bgq_xxcpnmadd(result_v2_c1, spinor1_v2_c1, qw, result_v2_c1);
	bgq_xmadd(result_v2_c2, qw, spinor1_v2_c2, spinor2_v2_c2);
	bgq_xxcpnmadd(result_v2_c2, spinor1_v2_c2, qw, result_v2_c2);

	bgq_xmadd(result_v3_c0, qw, spinor1_v3_c0, spinor2_v3_c0);
	bgq_xxcpnmadd(result_v3_c0, spinor1_v3_c0, qw, result_v3_c0);
	bgq_xmadd(result_v3_c1, qw, spinor1_v3_c1, spinor2_v3_c1);
	bgq_xxcpnmadd(result_v3_c1, spinor1_v3_c1, qw, result_v3_c1);
	bgq_xmadd(result_v3_c2, qw, spinor1_v3_c2, spinor2_v3_c2);
	bgq_xxcpnmadd(result_v3_c2, spinor1_v3_c2, qw, result_v3_c2);
#else
	bgq_su3_cjgvmul(imul1_v2, qc, spinor1_v2);
	bgq_su3_cjgvmul(imul1_v3, qc, spinor1_v3);
	bgq_su3_vsub(result_v2, spinor2_v2, imul1_v2);
	bgq_su3_vsub(result_v3, spinor2_v3, imul1_v3);
#endif

	bgq_su3_spinor_mov(*target, result);
}

#define OPERATOR_NAME bgq_spinorfield_cmul_plain_sub_gamma5_raw
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_cmul_plain_sub_gamma5
#define OPERATOR_EXTRAPARMS bgq_params(qc),bgq_params(qw)
#define OPERATOR_EXTRAARGS bgq_vars(qc),bgq_vars(qw)
#include "bgq_operator.inc.c"

void bgq_spinorfield_cmul_plain_sub_gamma5_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield1, bgq_weylfield_controlblock *sourcefield2, complex_double c) {
	bgq_vector4double_decl(qc);
	bgq_cconst(qc, creal(c), cimag(c));
	bgq_vector4double_decl(qw);
	bgq_neg(qw, qc);
	bgq_spinorfield_cmul_plain_sub_gamma5_raw_double(targetfield, isOdd, sourcefield1, sourcefield2, bgq_vars(qc), bgq_vars(qw));
}

void bgq_spinorfield_cmul_plain_sub_gamma5_float(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield1, bgq_weylfield_controlblock *sourcefield2, complex_double c) {
	bgq_vector4double_decl(qc);
	bgq_cconst(qc, creal(c), cimag(c));
	bgq_vector4double_decl(qw);
	bgq_neg(qw, qc);
	bgq_spinorfield_cmul_plain_sub_gamma5_raw_float(targetfield, isOdd, sourcefield1, sourcefield2, bgq_vars(qc), bgq_vars(qw));
}



static inline void bgq_site_cjgmul(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor), bgq_params(qc), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_su3_cvmul(result_v0,qc,spinor_v0);
	bgq_su3_cvmul(result_v1,qc,spinor_v1);
	bgq_su3_cjgvmul(result_v2,qc,spinor_v2);
	bgq_su3_cjgvmul(result_v3,qc,spinor_v3);

	bgq_su3_spinor_mov(*target, result);
}

#define OPERATOR_NAME bgq_spinorfield_cjgmul_raw
#define OPERATOR_ARGFIELDS 1
#define OPERATOR_VECSITEFUNC bgq_site_cjgmul
#define OPERATOR_EXTRAPARMS bgq_params(qc)
#define OPERATOR_EXTRAARGS bgq_vars(qc)
#include "bgq_operator.inc.c"

void bgq_spinorfield_cjgmul_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield, complex_double c) {
	bgq_vector4double_decl(qc);
	bgq_cconst(qc, creal(c), cimag(c));
	bgq_spinorfield_cjgmul_raw_double(targetfield, isOdd, sourcefield, bgq_vars(qc));
}

void bgq_spinorfield_cjgmul_float(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield, complex_double c) {
	bgq_vector4double_decl(qc);
	bgq_cconst(qc, creal(c), cimag(c));
	bgq_spinorfield_cjgmul_raw_float(targetfield, isOdd, sourcefield, bgq_vars(qc));
}



static inline void bgq_site_rmul_rmul_add(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), bgq_params(qr1), bgq_params(qr2), ucoord ic) {
	// The compiler will join this to FMA-instructions
	bgq_su3_spinor_decl(rmul1);
	bgq_su3_spinor_mul(rmul1, qr1, spinor1);

	bgq_su3_spinor_decl(rmul2);
	bgq_su3_spinor_mul(rmul2, qr2, spinor2);

	bgq_su3_spinor_decl(result);
	bgq_su3_spinor_add(result, rmul1, rmul2);

	bgq_su3_spinor_mov(*target, result);
}

#define OPERATOR_NAME bgq_spinorfield_rmul_rmul_add_raw
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_rmul_rmul_add
#define OPERATOR_EXTRAPARMS bgq_params(qr1), bgq_params(qr2)
#define OPERATOR_EXTRAARGS bgq_vars(qr1), bgq_vars(qr2)
#include "bgq_operator.inc.c"


void bgq_spinorfield_rmul_rmul_add_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield1, bgq_weylfield_controlblock *sourcefield2, double r1, double r2) {
	bgq_vector4double_decl(qr1);
	bgq_cconst(qr1, r1, r1);
	bgq_vector4double_decl(qr2);
	bgq_cconst(qr2, r2, r2);
	bgq_spinorfield_rmul_rmul_add_raw_double(targetfield, isOdd, sourcefield1, sourcefield2, bgq_vars(qr1), bgq_vars(qr2));
}

void bgq_spinorfield_rmul_rmul_add_float(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield1, bgq_weylfield_controlblock *sourcefield2, double r1, double r2) {
	bgq_vector4double_decl(qr1);
	bgq_cconst(qr1, r1, r1);
	bgq_vector4double_decl(qr2);
	bgq_cconst(qr2, r2, r2);
	bgq_spinorfield_rmul_rmul_add_raw_float(targetfield, isOdd, sourcefield1, sourcefield2, bgq_vars(qr1), bgq_vars(qr2));
}





static inline void bgq_site_cmul_plain_add(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), bgq_params(qc), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_su3_cvmadd(result_v0, qc, spinor1_v0, spinor2_v0);
	bgq_su3_cvmadd(result_v1, qc, spinor1_v1, spinor2_v1);
	bgq_su3_cvmadd(result_v2, qc, spinor1_v2, spinor2_v2);
	bgq_su3_cvmadd(result_v3, qc, spinor1_v3, spinor2_v3);

	bgq_su3_spinor_mov(*target, result);
}

#define OPERATOR_NAME bgq_spinorfield_cmul_plain_add_raw
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_cmul_plain_add
#define OPERATOR_EXTRAPARMS bgq_params(qc)
#define OPERATOR_EXTRAARGS bgq_vars(qc)
#include "bgq_operator.inc.c"

void bgq_spinorfield_cmul_plain_add_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield1, bgq_weylfield_controlblock *sourcefield2, complex_double c) {
	bgq_vector4double_decl(qc);
	bgq_cconst(qc, creal(c), cimag(c));
	bgq_spinorfield_cmul_plain_add_raw_double(targetfield, isOdd, sourcefield1, sourcefield2, bgq_vars(qc));
}

void bgq_spinorfield_cmul_plain_add_float(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield1, bgq_weylfield_controlblock *sourcefield2, complex_double c) {
	bgq_vector4double_decl(qc);
	bgq_cconst(qc, creal(c), cimag(c));
	bgq_spinorfield_cmul_plain_add_raw_float(targetfield, isOdd, sourcefield1, sourcefield2, bgq_vars(qc));
}
