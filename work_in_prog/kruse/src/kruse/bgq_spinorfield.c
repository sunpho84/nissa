/*
 * bgq_spinorfield.c
 *
 *  Created on: Oct 26, 2012
 *      Author: meinersbur
 */

#define BGQ_SPINORFIELD_C_
#include "bgq_spinorfield.h"

#include "bgq_HoppingMatrix.h"
#include "bgq_qpx.h"
#include "bgq_dispatch.h"
#include "bgq_comm.h"
#include "bgq_workers.h"

#include "../interface/internal_interface.h"

//#include "../geometry_eo.h"
//#include "../read_input.h"

#include <execinfo.h>
#include <mpi.h>
#include <sys/stat.h>
#include <stddef.h>
#ifndef NVALGRIND
#include <valgrind/memcheck.h>
#endif
#include <math.h>


typedef struct {
	double _Complex c[3];
} su3_vector_array64;

typedef struct {
	su3_vector_array64 v[4];
} spinor_array64;


double bgq_spinorfield_compare(bool isOdd, bgq_weylfield_controlblock *bgqfield, bgq_weylfield_controlblock *reffield, bool silent) {
	assert(bgqfield);
	assert(reffield);
	//assert(bgqfield->isInitialized);
	assert(bgqfield->isOdd == isOdd);

	//bool readFulllayout = bgqfield->hasFullspinorData;
	//bgq_spinorfield_layout layout = bgq_spinorfield_prepareRead(bgqfield, isOdd, true, true, true, true);
	//bgq_spinorfield_setup(bgqfield, isOdd, readFulllayout, false, !readFulllayout, false, false);
	//bgq_master_sync(); // Necessary after bgq_spinorfield_setup if field is accessed without bgq_master_call (which does this implicitely)

	double diff_max = 0;
	size_t count = 0;
	for (size_t z = 0; z < LOCAL_LZ ; z += 1) {
		for (size_t y = 0; y < LOCAL_LY ; y += 1) {
			for (size_t x = 0; x < LOCAL_LX ; x += 1) {
				for (size_t t = 0; t < LOCAL_LT ; t += 1) {
					if (bgq_local2isOdd(t, x, y, z) != isOdd)
						continue;

					bgq_spinor bgqspinor = bgq_spinorfield_getspinor(bgqfield, t,x,y,z);
					bgq_spinor refspinor = bgq_spinorfield_getspinor(reffield, t,x,y,z);

					bool first = true;
					for (size_t v = 0; v < 4; v += 1) {
						for (size_t c = 0; c < 3; c += 1) {
							complexdouble bgqvalue = bgqspinor.v[v].c[c];
							complexdouble refvalue = refspinor.v[v].c[c];

							double diff = cabs(bgqvalue - refvalue);

							if (diff > 0.01 || isnan(diff)) {
								if (!silent && first) {
									master_print("Coordinate (%zu,%zu,%zu,%zu)(%zu,%zu): ref=(%f + %fi) != bgq=(%f + %fi) off by %f\n", t, x, y, z, v, c, creal(refvalue), cimag(refvalue), creal(bgqvalue), cimag(bgqvalue), diff);
								}
								if (first)
									count += 1;
								first = false;
							}
							if (diff > diff_max || isnan(diff)) {
								diff_max = diff;
							}
						}
					}
				}
			}
		}
	}

	double global_diff_max;
	MPI_Allreduce(&diff_max, &global_diff_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	if (count > 0) {
		if (!silent)
			master_print(" %zu sites of %d wrong\n", count, VOLUME/2);
	}

	return global_diff_max;
}


static double bgq_spinorfield_legacy_compare(bool isOdd, bgq_weylfield_controlblock *bgqfield, bgq_spinorfield_layout bgqlayout, spinor *reffield, bool silent) {
	assert(bgqfield);
	assert(reffield);
	//assert(bgqfield->isInitialized);
	assert(bgqfield->isOdd == isOdd);

	//bool readFulllayout = bgqfield->hasFullspinorData;
	//bgq_spinorfield_layout layout = bgq_spinorfield_prepareRead(bgqfield, isOdd, true, true, true, true);
	//bgq_spinorfield_setup(bgqfield, isOdd, readFulllayout, false, !readFulllayout, false, false);
	bgq_master_sync(); // Necessary after bgq_spinorfield_setup if field is accessed without bgq_master_call (which does this implicitely)

	double diff_max = 0;
	size_t count = 0;
	for (size_t z = 0; z < LOCAL_LZ ; z += 1) {
		for (size_t y = 0; y < LOCAL_LY ; y += 1) {
			for (size_t x = 0; x < LOCAL_LX ; x += 1) {
				for (size_t t = 0; t < LOCAL_LT ; t += 1) {
					if (bgq_local2isOdd(t, x, y, z) != isOdd)
						continue;

					const int ix = Index(t,x,y,z); /* lexic coordinate */
					assert(ix == Index(t,x,y,z));
					//int iy = g_lexic2eo[ix]; /* even/odd coordinate (even and odd sites in two different fields of size VOLUME/2, first even field followed by odd) */
					//assert(0 <= iy && iy < (VOLUME+RAND));
					int icx = g_lexic2eosub[ix]; /*  even/odd coordinate relative to field base */
					assert(0 <= icx && icx < VOLUME/2);
					//assert(icx == iy - (isOdd ? (VOLUME+RAND)/2 : 0));
					spinor_array64 *sp = (spinor_array64*) &reffield[icx];

					size_t ic = bgq_local2collapsed(t, x, y, z);
					size_t k = bgq_local2k(t, x, y, z);
					//bgq_spinor bgqspinor = bgq_spinorfield_getspinor(bgqfield, t,x,y,z);
					bgq_su3_spinor_decl(qpxspinor);
					bgq_spinorfield_readSpinor(&qpxspinor, bgqfield, isOdd, ic, bgqlayout&ly_weyl, bgqlayout&ly_sloppy, bgqlayout&ly_mul, bgqlayout==ly_legacy);
					bgq_spinor bgqspinor = bgq_spinor_fromqpx(qpxspinor, k);

					bool first = true;
					for (size_t v = 0; v < 4; v += 1) {
						for (size_t c = 0; c < 3; c += 1) {
							complexdouble bgqvalue = bgqspinor.v[v].c[c];
							if (isnan(creal(bgqvalue)) || isnan(cimag(bgqvalue))) {
								int b = 0;
							}
							complexdouble refvalue = sp->v[v].c[c];
							if (isnan(creal(refvalue)) || isnan(cimag(refvalue))) {
								int a = 0;
							}

							double diff = cabs(bgqvalue - refvalue);

							if (diff > 0.01 || isnan(diff)) {
								if (!silent && first)
									master_print("Coordinate (%zu,%zu,%zu,%zu)(%zu,%zu): ref=(%f + %fi) != bgq=(%f + %fi) off by %f\n", t, x, y, z, v, c, creal(refvalue), cimag(refvalue), creal(bgqvalue), cimag(bgqvalue), diff);
								if (first)
									count += 1;
								first = false;
							}
 							if (diff > diff_max || isnan(diff)) {
								diff_max = diff;
							}
						}
					}
				}
			}
		}
	}

	double global_diff_max;
	MPI_Allreduce(&diff_max, &global_diff_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	if (count > 0) {
		if (!silent)
			master_print(" %zu sites of %d wrong\n", count, VOLUME/2);
	}

	return global_diff_max;
}


void bgq_weyl_expect(bgq_weyl_nonvec weyl, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
	size_t t_global = bgq_local2global_t(t);
	size_t x_global = bgq_local2global_x(x);
	size_t y_global = bgq_local2global_y(y);
	size_t z_global = bgq_local2global_z(z);
	if (isSrc) {
		bgq_direction_move_global(&t_global, &x_global, &y_global, &z_global, d);
		d = bgq_direction_revert(d);
	}

	assert(weyl.s[0][0] == t_global);
	assert(weyl.s[0][1] == x_global);
	assert(weyl.s[0][2] == y_global);
	assert(weyl.s[1][0] == z_global);
	assert(weyl.s[1][1] == d);
	assert(weyl.s[1][2] == 0.125);
}


static bgq_weyl_vec_double *bgq_offset2pointer_double(bgq_weylfield_controlblock *field, size_t offset) {
	assert(field);
	assert(offset % sizeof(bgq_weyl_vec_double) == 0);
	assert(bgq_weyl_section_offset(0) <= offset && offset < bgq_weyl_section_offset(sec_end));

	bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
	bgq_weyl_vec_double *baseptr = bgq_section_baseptr_double(field, sec);
	assert(baseptr);
	size_t baseoffset = bgq_weyl_section_offset(sec);
	assert(offset >= baseoffset);
	size_t reloffset = offset - baseoffset;
	return (bgq_weyl_vec_double*)((uint8_t*)baseptr + reloffset);
}


static bgq_weyl_vec_float *bgq_offset2pointer_float(bgq_weylfield_controlblock *field, size_t offset) {
	assert(field);
	assert(offset % sizeof(bgq_weyl_vec_double) == 0);
	assert(bgq_weyl_section_offset(0) <= offset && offset < bgq_weyl_section_offset(sec_end));

	bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
	bgq_weyl_vec_float *baseptr = bgq_section_baseptr_float(field, sec);
	assert(baseptr);
	size_t baseoffset = bgq_weyl_section_offset(sec);
	assert(offset >= baseoffset);
	size_t reloffset = offset - baseoffset;
	return (bgq_weyl_vec_float*)((uint8_t*)baseptr + reloffset/2);
}


static bgq_weyl_vec_double *bgq_index2pointer_double(bgq_weylfield_controlblock *field, ucoord index) {
	return bgq_offset2pointer_double(field, bgq_index2offset(index));
}

static bgq_weyl_vec_float *bgq_index2pointer_float(bgq_weylfield_controlblock *field, ucoord index) {
	return bgq_offset2pointer_float(field, bgq_index2offset(index));
}


static size_t bgq_weylfield_bufferoffset2consecutiveoffset(bool isOdd, size_t offset, size_t k) {
	assert(offset);
	assert(offset % sizeof(bgq_weyl_vec_double) == 0);
	bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
	size_t index = bgq_offset2index(offset);

	//bgq_direction d = bgq_section2direction(sec);
	//assert(d == bgq_offset2ddst(offset));
	bgq_direction d = bgq_offset2ddst(offset);
	ucoord ic = g_bgq_index2collapsed[isOdd][index];
	assert(ic != -1);
	return bgq_collapsed2consecutiveoffset(ic, d);
}


static bgq_spinor_nonvec bgq_spinor_coord_encode(scoord t, scoord x, scoord y, scoord z) {
	ucoord t_global = bgq_local2global_t(t);
	ucoord x_global = bgq_local2global_x(x);
	ucoord y_global = bgq_local2global_y(y);
	ucoord z_global = bgq_local2global_z(z);

	bgq_spinor_nonvec result = {{{{0}}}};
	result.v[0].c[0] = t_global;
	result.v[0].c[1] = x_global;
	result.v[0].c[2] = y_global;
	result.v[1].c[0] = z_global;

	return result;
}


static void bgq_spinorveck_write_double(bgq_spinorsite_double *target, ucoord k, bgq_spinor data) {
	for (ucoord i = 0; i < 4; i+=1) {
		for (ucoord l = 0; l < 3; l+=1) {
			target->s[i][l][k] = data.v[i].c[l];
		}
	}
}


static void bgq_spinorveck_write_float(bgq_spinorsite_float *target, ucoord k, bgq_spinor data) {
	for (ucoord i = 0; i < 4; i+=1) {
		for (ucoord l = 0; l < 3; l+=1) {
			target->s[i][l][k] = data.v[i].c[l];
		}
	}
}


void bgq_spinorveck_written_double(bgq_spinorsite_double *targetspinor, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z) {
	assert(targetspinor->s[0][0][1] == targetspinor->s[0][0][1]); // For valgrind (or to catch NaN)
#ifdef BGQ_COORDCHECK
	bgq_spinor coord = bgq_spinor_coord_encode(t,x,y,z);
	bgq_spinorveck_write_double(targetspinor, k, coord);
#endif
}

void bgq_spinorveck_written_float(bgq_spinorsite_float *targetspinor, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z) {
	if (targetspinor->s[0][0][0] == 0)
		assert(targetspinor->s[0][0][1] != 1); // For valgrind
#ifdef BGQ_COORDCHECK
	printf("ciao\n");
	bgq_spinor coord = bgq_spinor_coord_encode(t,x,y,z);
	bgq_spinorveck_write_float(targetspinor, k, coord);
#endif
}

static bgq_weyl_nonvec bgq_weyl_coord_encode(ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
	ucoord t_global = bgq_local2global_t(t);
	ucoord x_global = bgq_local2global_x(x);
	ucoord y_global = bgq_local2global_y(y);
	ucoord z_global = bgq_local2global_z(z);
	if (isSrc) {
		bgq_direction_move_global(&t_global, &x_global, &y_global, &z_global, d);
		d = bgq_direction_revert(d);
	}

	bgq_weyl_nonvec result = {{{0}}};
	result.s[0][0] = t_global;
	result.s[0][1] = x_global;
	result.s[0][2] = y_global;
	result.s[1][0] = z_global;
	result.s[1][1] = d;
	result.s[1][2] = 0.125;
	return result;
}

static void bgq_weylveck_write_double(bgq_weyl_vec_double *target, ucoord k, bgq_weyl_nonvec data) {
	for (ucoord i = 0; i < 2; i += 1) {
		for (ucoord l = 0; l < 3; l += 1) {
			target->s[i][l][k] = data.s[i][l];
		}
	}
}


static void bgq_weylveck_write_float(bgq_weyl_vec_float *target, ucoord k, bgq_weyl_nonvec data) {
	for (ucoord i = 0; i < 2; i += 1) {
		for (ucoord l = 0; l < 3; l += 1) {
			target->s[i][l][k] = data.s[i][l];
		}
	}
}

#define bgq_weylveck_written NAME2(bgq_weylveck_written,PRECISION)
void bgq_weylveck_written_double(bgq_weyl_vec_double *targetweyl, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
#ifdef BGQ_COORDCHECK
  //printf("writing tloc: %d, tglob: %d\n",t,bgq_local2global_t(t));
	bgq_weyl_nonvec coord = bgq_weyl_coord_encode(t,x,y,z,d,isSrc);
	bgq_weylveck_write_double(targetweyl, k, coord);
#endif
}


void bgq_weylveck_written_float(bgq_weyl_vec_float *targetweyl, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
#ifdef BGQ_COORDCHECK
	bgq_weyl_nonvec coord = bgq_weyl_coord_encode(t,x,y,z,d,isSrc);
	bgq_weylveck_write_float(targetweyl, k, coord);
#endif
}


static size_t bgq_physical_halo_sites(bgq_dimension dim) {
	switch (dim) {
	case DIM_T:
		return HALO_T * LOCAL_HALO_T/PHYSICAL_LP;
	case DIM_X:
		return PHYSICAL_HALO_X;
	case DIM_Y:
		return PHYSICAL_HALO_Y;
	case DIM_Z:
		return PHYSICAL_HALO_Z;
	}
	UNREACHABLE
	return -1;
}


static size_t bgq_physical_halo_vecs(bgq_dimension dim) {
	switch (dim) {
	case DIM_T:
		return HALO_T * LOCAL_HALO_T/(PHYSICAL_LP*PHYSICAL_LK);
	case DIM_X:
		return PHYSICAL_HALO_X;
	case DIM_Y:
		return PHYSICAL_HALO_Y;
	case DIM_Z:
		return PHYSICAL_HALO_Z;
	}
	UNREACHABLE
	return -1;
}


static bgq_weylfield_section bgq_direction2writesec(bgq_direction d) {
	switch (d) {
	case TUP:
		if (BGQ_UNVECTORIZE || !COMM_T) {
			return sec_temp_tup;
		} else {
			return sec_send_tup;
		}
	case TDOWN:
		if (BGQ_UNVECTORIZE || !COMM_T) {
			return sec_temp_tdown;
		} else {
			return sec_send_tdown;
		}
	case XUP:
		return sec_send_xup;
	case XDOWN:
		return sec_send_xdown;
	case YUP:
		return sec_send_yup;
	case YDOWN:
		return sec_send_ydown;
	case ZUP:
		return sec_send_zup;
	case ZDOWN:
		return sec_send_zdown;
	}
	UNREACHABLE
	return -1;
}


void bgq_spinorfield_enableLayout(bgq_weylfield_controlblock *field, tristate isOdd, bgq_spinorfield_layout layout, bool disableOthers, bool preserveData) {
	assert(field);
	assert(!field->pendingDatamove);

	if (preserveData) {
		 if (field->isOdd==tri_unknown) {
		 } else if (isOdd==tri_unknown) {
			isOdd = field->isOdd;
		} else {
			assert(isOdd==field->isOdd);
		}
	}

	bool isSloppy = layout&ly_sloppy;

	// Allocate necessary fields
	if (layout==ly_legacy) {
		// Allocation has been done by init_spinor_field()
		assert(field->legacy_field);
	} else if (layout & ly_weyl) {
		assert(isOdd!=tri_unknown);

		// == Fields itself ==
		if (!field->sec_collapsed_double) {
			field->sec_collapsed_double = malloc_aligned(PHYSICAL_VOLUME * sizeof(*field->sec_collapsed_double), BGQ_ALIGNMENT_L2);
#ifndef NVALGRIND
			VALGRIND_CREATE_MEMPOOL(field->sec_collapsed_double, 0, false);
			VALGRIND_MEMPOOL_ALLOC(field->sec_collapsed_double, field->sec_collapsed_double, PHYSICAL_VOLUME * sizeof(*field->sec_collapsed_double));
#endif
		}

		if (!field->sec_collapsed_float) {
			field->sec_collapsed_float = (bgq_weylsite_float*)field->sec_collapsed_double;
		}


		// == Pointer into fields ==
		if ((isSloppy && !field->sendptr_float[isOdd]) || (!isSloppy && !field->sendptr_double[isOdd])) {
			if (isSloppy) {
				field->sendptr_float[isOdd] = malloc_aligned(PHYSICAL_VOLUME * sizeof(*field->sendptr_float[isOdd]), BGQ_ALIGNMENT_L2);
			} else {
				field->sendptr_double[isOdd] = malloc_aligned(PHYSICAL_VOLUME * sizeof(*field->sendptr_double[isOdd]), BGQ_ALIGNMENT_L2);
			}
			for (size_t d = 0; d < PHYSICAL_LD; d+=1) {
				bgq_dimension dim = bgq_direction2dimension(d);

				if (isSloppy) {
					field->consptr_float[isOdd][d] = malloc_aligned(bgq_physical_halo_sites(dim) * sizeof(*field->consptr_float[isOdd][d]), BGQ_ALIGNMENT_L2);
				} else {
					field->consptr_double[isOdd][d] = malloc_aligned(bgq_physical_halo_sites(dim) * sizeof(*field->consptr_double[isOdd][d]), BGQ_ALIGNMENT_L2);
				}
			}


			// For main kernel (surface & body)
			for (ucoord ic_src = 0; ic_src < PHYSICAL_VOLUME; ic_src += 1) {
				for (size_t d_dst = 0; d_dst < PHYSICAL_LD; d_dst += 1) {
					ucoord index = g_bgq_collapsed2indexsend[isOdd/*_dst*/][ic_src].d[d_dst];

					if (isSloppy) {
						bgq_weyl_vec_float *ptr = bgq_index2pointer_float(field, index);
						field->sendptr_float[isOdd][ic_src].d[d_dst] = ptr;
					} else {
						bgq_weyl_vec_double *ptr = bgq_index2pointer_double(field, index);
						field->sendptr_double[isOdd][ic_src].d[d_dst] = ptr;
					}
				}
			}


			// For 5th phase (datamove)
			for (ucoord d_dst = 0; d_dst < PHYSICAL_LD; d_dst+=1) {
				bgq_direction d_src = bgq_direction_revert(d_dst);
				bgq_dimension dim = bgq_direction2dimension(d_dst);
				ucoord sites = bgq_physical_halo_sites(dim);
				for (ucoord j = 0; j < sites; j+=1) {
					bgq_weylfield_section sec = bgq_direction2writesec(d_src);
					size_t baseoffset = bgq_weyl_section_offset(sec);
					ucoord baseindex = bgq_offset2index(baseoffset);
					ucoord index = baseindex + j;
					ucoord ic_dst = g_bgq_index2collapsed[isOdd][index]; // Found out what the previous phase wrote here
					size_t offset_cons = bgq_collapsed2consecutiveoffset(ic_dst, d_dst);

					if (isSloppy) {
						bgq_weyl_vec_float *ptr = bgq_offset2pointer_float(field, offset_cons);
						field->consptr_float[isOdd][d_dst][j] = ptr;
					} else {
						bgq_weyl_vec_double *ptr = bgq_offset2pointer_double(field, offset_cons);
						field->consptr_double[isOdd][d_dst][j] = ptr;
					}
				}
			}
		}
	} else {
		if (!field->sec_fullspinor_double) {
			field->sec_fullspinor_double = malloc_aligned(LOCAL_VOLUME/PHYSICAL_LP * sizeof(*field->sec_fullspinor_double), BGQ_ALIGNMENT_L2);
#ifndef NVALGRIND
			VALGRIND_CREATE_MEMPOOL(field->sec_fullspinor_double, 0, false);
			VALGRIND_MEMPOOL_ALLOC(field->sec_fullspinor_double, field->sec_fullspinor_double, LOCAL_VOLUME/PHYSICAL_LP * sizeof(*field->sec_fullspinor_double));
#endif
		}
		if (!field->sec_fullspinor_float) {
			field->sec_fullspinor_float = (bgq_spinorsite_float*)field->sec_fullspinor_double;
		}
	}


	if (disableOthers) {
		field->has_fulllayout_double = false;
		field->has_fulllayout_float = false;
		field->has_weyllayout_double = false;
		field->has_weyllayout_float = false;
		field->has_legacy = false;
	}
	field->isOdd = isOdd;


	switch (layout) {
	case ly_full_double:
		assert(field->sec_fullspinor_double);
		if (!preserveData) {
#ifndef NDEBUG
		memset(field->sec_fullspinor_double, 0xFF, LOCAL_VOLUME/PHYSICAL_LP * sizeof(*field->sec_fullspinor_double));
#endif
#ifndef NVALGRIND
		VALGRIND_MEMPOOL_FREE(field->sec_fullspinor_double, field->sec_fullspinor_double);
		VALGRIND_MEMPOOL_ALLOC(field->sec_fullspinor_double, field->sec_fullspinor_double, LOCAL_VOLUME/PHYSICAL_LP * sizeof(*field->sec_fullspinor_double));
#endif
		}
		field->has_fulllayout_double = true;
		break;
	case ly_full_float:
		assert(field->sec_fullspinor_float);
		if (!preserveData) {
#ifndef NDEBUG
		memset(field->sec_fullspinor_float, 0xFF, LOCAL_VOLUME/PHYSICAL_LP * sizeof(*field->sec_fullspinor_float));
#endif
#ifndef NVALGRIND
		VALGRIND_MEMPOOL_FREE(field->sec_fullspinor_float, field->sec_fullspinor_float);
		VALGRIND_MEMPOOL_ALLOC(field->sec_fullspinor_float, field->sec_fullspinor_float, LOCAL_VOLUME/PHYSICAL_LP * sizeof(*field->sec_fullspinor_float));
#endif
		}
		field->has_fulllayout_float = true;
		break;
	case ly_weyl_double:
		assert(field->sec_collapsed_double);
		if (!preserveData) {
#ifndef NDEBUG
		memset(field->sec_collapsed_double, 0xFF, PHYSICAL_VOLUME * sizeof(*field->sec_collapsed_double));
#endif
#ifndef NVALGRIND
		VALGRIND_MEMPOOL_FREE(field->sec_collapsed_double, field->sec_collapsed_double);
		VALGRIND_MEMPOOL_ALLOC(field->sec_collapsed_double, field->sec_collapsed_double, PHYSICAL_VOLUME * sizeof(*field->sec_collapsed_double));
#endif
		}
		field->has_weyllayout_double = true;
		break;
	case ly_weyl_float:
		assert(field->sec_collapsed_float);
		if (!preserveData) {
#ifndef NDEBUG
		memset(field->sec_collapsed_float, 0xFF, PHYSICAL_VOLUME * sizeof(*field->sec_collapsed_float));
#endif
#ifndef NVALGRIND
		VALGRIND_MEMPOOL_FREE(field->sec_collapsed_float, field->sec_collapsed_float);
		VALGRIND_MEMPOOL_ALLOC(field->sec_collapsed_float, field->sec_collapsed_float, PHYSICAL_VOLUME * sizeof(*field->sec_collapsed_float));
#endif
		}
		field->has_weyllayout_float = true;
		break;
	case ly_legacy:
		assert(field->legacy_field);
		if (!preserveData) {
#ifndef NDEBUG
		memset(field->legacy_field, 0xFF, VOLUMEPLUSRAND/2 * sizeof(*field->legacy_field));
#endif
#ifndef NVALGRIND
		VALGRIND_MEMPOOL_FREE(field->collectionBase->legacy_base, field->legacy_field);
		VALGRIND_MEMPOOL_ALLOC(field->collectionBase->legacy_base, field->legacy_field, VOLUMEPLUSRAND/2*sizeof(*field->legacy_field));
#endif
		}
		field->has_legacy = true;
		break;
	default:
		assert(!"Not yet implemented");
		break;
	}
}


bgq_spinor bgq_legacy_getspinor(spinor *spinor, ucoord t, ucoord x, ucoord y, ucoord z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	spinorfield_enable(spinor, true, false);
	int ix = Index(t,x,y,z); /* lexic coordinate */
	assert(ix == Index(t,x,y,z));
	//int iy = g_lexic2eo[ix]; /* even/odd coordinate (even and odd sites in two different fields of size VOLUME/2, first even field followed by odd) */
	//assert(0 <= iy && iy < (VOLUME+RAND));
	int icx = g_lexic2eosub[ix]; /*  even/odd coordinate relative to field base */
	assert(0 <= icx && icx < VOLUME/2);

	bgq_spinor *result = (bgq_spinor*)&spinor[icx];
	return *result;
}


bgq_spinor bgq_spinorfield_getspinor(bgq_weylfield_controlblock *field, ucoord t, ucoord x, ucoord y, ucoord z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);

	//assert(bgq_local2isOdd(t, x, y, z) == field->isOdd);
	bool isOdd = bgq_local2isOdd(t, x, y, z);
	ucoord ic = bgq_local2collapsed(t, x, y, z);
	ucoord k = bgq_local2k(t, x, y, z);

	bgq_spinorfield_layout layout = bgq_spinorfield_prepareRead(field, isOdd, true, true, true, true, true);
	//bgq_spinorfield_setup(field, field->isOdd, !(layout & ly_weyl), false, (layout & ly_weyl), false, false);
	bgq_master_sync();
	bgq_su3_spinor_decl(spinor);
	bgq_spinorfield_readSpinor(&spinor, field, isOdd, ic, layout&ly_weyl, layout&ly_sloppy, layout&ly_mul, layout==ly_legacy);
	return bgq_spinor_fromqpx(spinor, k);
}


char *(g_idxdesc[BGQREF_count]);
complexdouble *g_bgqvalue = NULL;
complexdouble *g_refvalue = NULL;
bool *g_bgqhasvalue = NULL;
bool *g_refhasvalue = NULL;

void bgq_initbgqref_impl(void) {
	size_t datasize = sizeof(complexdouble) * VOLUME * lengthof(g_idxdesc);
	if (g_refvalue == NULL) {
		g_bgqvalue = malloc_aligned(datasize, BGQ_ALIGNMENT_L2);
		g_refvalue = malloc_aligned(datasize, BGQ_ALIGNMENT_L2);
		g_bgqhasvalue = malloc(VOLUME * lengthof(g_idxdesc) * sizeof(*g_bgqhasvalue));
		g_refhasvalue = malloc(VOLUME * lengthof(g_idxdesc) * sizeof(*g_refhasvalue));
	}
	memset(g_bgqvalue, 0xFF, datasize);
	memset(g_refvalue, 0xFF, datasize);
	memset(g_refhasvalue, 0, VOLUME * lengthof(g_idxdesc) * sizeof(*g_refhasvalue));
	memset(g_bgqhasvalue, 0, VOLUME * lengthof(g_idxdesc) * sizeof(*g_bgqhasvalue));

	for (int idx = 0; idx <  lengthof(g_idxdesc); idx+=1) {
		g_idxdesc[idx] = NULL;
	}
}


void bgq_setdesc_impl(int idx, char *desc){
	g_idxdesc[idx] = desc;
}


void bgq_setrefvalue_impl(int t, int x, int y, int z, bgqref idx, complexdouble val) {
	if (t < 0)
		t = 0;
	if (t >= LOCAL_LT)
		t = LOCAL_LT-1;
	if (x < 0)
		x = 0;
	if (x >= LOCAL_LX)
		x = LOCAL_LX-1;
	if (y < 0)
		y = 0;
	if (y >= LOCAL_LY)
		y = LOCAL_LY-1;
	if (z < 0)
		z = 0;
	if (z >= LOCAL_LZ)
		z = LOCAL_LZ-1;
	g_refvalue[(((idx*LOCAL_LT + t)*LOCAL_LX + x)*LOCAL_LY + y)*LOCAL_LZ + z] = val;
	g_refhasvalue[(((idx*LOCAL_LT + t)*LOCAL_LX + x)*LOCAL_LY + y)*LOCAL_LZ + z] = true;
	if (!g_idxdesc[idx])
		g_idxdesc[idx] = "";
}


void bgq_setbgqvalue_impl(int t, int x, int y, int z, bgqref idx, complex_double val) {
	if (idx==BGQREF_TUP && t==0 && x==0 && y==0 && z==0) {
		int a = 0;
	}
	if (val == 0) {
		if (val == 1) {
			//assert(false);
		}
	}

	if (t < 0)
		t = 0;
	if (t >= LOCAL_LT)
		t = LOCAL_LT-1;
	if (x < 0)
		x = 0;
	if (x >= LOCAL_LX)
		x = LOCAL_LX-1;
	if (y < 0)
		y = 0;
	if (y >= LOCAL_LY)
		y = LOCAL_LY-1;
	if (z < 0)
		z = 0;
	if (z >= LOCAL_LZ)
		z = LOCAL_LZ-1;
	g_bgqvalue[(((idx*LOCAL_LT + t)*LOCAL_LX + x)*LOCAL_LY + y)*LOCAL_LZ + z] = val;
	g_bgqhasvalue[(((idx*LOCAL_LT + t)*LOCAL_LX + x)*LOCAL_LY + y)*LOCAL_LZ + z] = true;
	if (!g_idxdesc[idx])
		g_idxdesc[idx] = "";
}


void bgq_setbgqvalue_src_impl(ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bgqref idx, complexdouble val) {
	bgq_direction_move_local(&t, &x, &y, &z, d);
	bgq_setbgqvalue(t, x, y, z, idx, val);
}


void bgq_savebgqref_impl(void) {
	if (g_proc_id != 0)
		return;

	int i = 0;
	while (true) {
		char filename[100];
		snprintf(filename, sizeof(filename)-1, "cmp_%d.txt", i);

		struct stat buf;
		if (stat(filename, &buf) != -1) {
			master_print("MK file %s already exists\n", filename);
			i += 1;
			continue;
		}

		// Create marker file
		FILE *file =  fopen(filename, "w");
		fclose(file);

		break;
	}

	for (int idx = 0; idx <  lengthof(g_idxdesc); idx+=1) {
		if (!g_idxdesc[idx])
			continue;

		char reffilename[100];
		char refdir[100];
		snprintf(refdir, sizeof(refdir)-1, "cmpref");
		snprintf(reffilename, sizeof(reffilename)-1, "%s/cmp_idx%02d_%s.txt", refdir, idx, g_idxdesc[idx]);

		char bgqdir[100];
		char bgqfilename[100];
		snprintf(bgqdir, sizeof(bgqdir)-1, "cmpbgq");
		snprintf(bgqfilename, sizeof(bgqfilename)-1, "%s/cmp_idx%02d_%s.txt", bgqdir, idx, g_idxdesc[idx]);
		master_print("Cmp Going to write to %s and %s\n", reffilename, bgqfilename);

		mkdir(refdir, S_IRWXU | S_IRWXG | S_IRWXO);
		mkdir(bgqdir, S_IRWXU | S_IRWXG | S_IRWXO);

		FILE *reffile = fopen(reffilename, "w");
		FILE *bgqfile = fopen(bgqfilename, "w");

		fprintf(reffile, "%s\n\n", g_idxdesc[idx]);
		fprintf(bgqfile, "%s\n\n", g_idxdesc[idx]);

		for (int z = 0; z < LOCAL_LZ; z += 1) {
			for (int x = 0; x < LOCAL_LX; x += 1) {
				for (int y = 0; y < LOCAL_LY; y += 1) {
					fprintf(reffile, "x=%d y=%d z=%d: ", x,y,z);
					fprintf(bgqfile, "x=%d y=%d z=%d: ", x,y,z);
					for (int t = 0; t < LOCAL_LT; t += 1) {
						complexdouble refval = g_refvalue[(((idx*LOCAL_LT + t)*LOCAL_LX + x)*LOCAL_LY + y)*LOCAL_LZ + z];
						complexdouble bgqval = g_bgqvalue[(((idx*LOCAL_LT + t)*LOCAL_LX + x)*LOCAL_LY + y)*LOCAL_LZ + z];
						bool refhasval = g_refhasvalue[(((idx*LOCAL_LT + t)*LOCAL_LX + x)*LOCAL_LY + y)*LOCAL_LZ + z];
						bool bgqhasval = g_bgqhasvalue[(((idx*LOCAL_LT + t)*LOCAL_LX + x)*LOCAL_LY + y)*LOCAL_LZ + z];

						if (refhasval) {
							fprintf(reffile, "%8f + %8fi	", creal(refval), cimag(refval));
						} else {
							fprintf(reffile, "                    	");
						}

						if (bgqhasval) {
							fprintf(bgqfile, "%8f + %8fi	", creal(bgqval), cimag(bgqval));
						} else {
							fprintf(bgqfile, "                    	");
						}
					}
					fprintf(reffile, "\n");
					fprintf(bgqfile, "\n");
				}
			}
		}

		fclose(reffile);
		fclose(bgqfile);

		master_print("Cmp data written to %s and %s\n", reffilename, bgqfilename);
	}

	bgq_initbgqref();
}


static void print_trace(FILE *target) {
	/* GG */
	if ( g_proc_id )
	return;

	void *array[50];
	size_t size;
	char **strings;
	size_t i;

	size = backtrace (array, 50);
	strings = backtrace_symbols (array, size);

	//printf ("Obtained %zd stack frames.\n", size);

	for (i = 0; i < size; i++) {
	fprintf (target, "%s\n", strings[i]);
	}

	free (strings);
}


void bgq_spinorfield_dump(bgq_weylfield_controlblock *field, char *desc) {
	assert(field);

	if (g_proc_id!=0)
		return;
	if (omp_get_thread_num()!=0)
		return;

	
	static int g_bgq_dumpnr = 0;
	int dumpnr = g_bgq_dumpnr;
	g_bgq_dumpnr += 1;

	char descSurrogate[100];
	if (!desc) {
		desc = "dump";
	}


	static char dir[100] = "";
	if (dir[0] == '\0') {
		int i = 0;
		while (true) {
			snprintf(dir, sizeof(dir)-1, "dump_%d", i);

			struct stat buf;
			if (stat(dir, &buf) != -1) {
				i += 1;
				continue;
			}

			mkdir(dir, S_IRWXU | S_IRWXG | S_IRWXO);
			break;
		}
	}
	
		char filename[100];
		int j = 0;
		while (true) {
			snprintf(filename, sizeof(filename)-1, "%s/spinorfield_%d_%d_%s.txt", dir, dumpnr, j, desc);

			struct stat buf;
			if (stat(filename, &buf) != -1) {
				j += 1;
				continue;
			}

			break;
		}

		master_print("Dumping to file %s\n", filename);
		//mkdir(dir, S_IRWXU | S_IRWXG | S_IRWXO);
		FILE *dfile = fopen(filename, "w");
		fprintf(dfile, "Sequence number %d\n", dumpnr);
		print_trace(dfile);
		fprintf(dfile, "\n\n");

		for (int z = 0; z < LOCAL_LZ; z += 1) {
			for (int x = 0; x < LOCAL_LX; x += 1) {
				for (int y = 0; y < LOCAL_LY; y += 1) {
					fprintf(dfile, "x=%d y=%d z=%d: ", x,y,z);
					for (int t = 0; t < LOCAL_LT; t += 1) {
						if (bgq_local2isOdd(t, x, y, z) == field->isOdd) {
							complexdouble val = bgq_spinorfield_getspinor(field, t, x, y, z).v[0].c[0];
							fprintf(dfile, "%8f + %8fi	", creal(val), cimag(val));
						} else {
							fprintf(dfile, "                    	");
						}
					}
					fprintf(dfile, "\n");
				}
			}
		}

		fclose(dfile);
		//master_print("Cmp data written to %s and %s\n", reffilename, bgqfilename);
}

void spinorfield_dump(const spinor *field, char *desc) {
	bgq_spinorfield_dump(bgq_translate_spinorfield(field), desc);
}


// not reliable!!!
static size_t bgq_fieldpointer2offset(void *ptr) {
	for (size_t i = 0; i < g_bgq_spinorfields_count; i+=1) {
		bgq_weylfield_controlblock *field = &g_bgq_spinorfields[i];
		size_t result = bgq_pointer2offset_raw(field, ptr, false);
		if (result != -1)
			return result;
	}

	assert(!"Pointer does not point to weyllayout");
	return -1;
}


bgq_weyl_vec_double *bgq_section_baseptr_double(bgq_weylfield_controlblock *field, bgq_weylfield_section section) {
	assert(field);
	//assert(field->isInitialized);
	bgq_weyl_vec_double *result = NULL;

	switch (section) {
	case sec_surface:
		if (PHYSICAL_SURFACE==0)
			return NULL;
		result = (bgq_weyl_vec_double*)&field->sec_collapsed_double[bgq_surface2collapsed(0)];
		break;
	case sec_body:
		if (PHYSICAL_BODY==0)
			return NULL;
		result = (bgq_weyl_vec_double*)&field->sec_collapsed_double[bgq_body2collapsed(0)];
		break;
	case sec_send_tup:
		result = g_bgq_sec_send_double[TUP];
		break;
	case sec_send_tdown:
		result = g_bgq_sec_send_double[TDOWN];
		break;
	case sec_send_xup:
		result = g_bgq_sec_send_double[XUP];
		break;
	case sec_send_xdown:
		result = g_bgq_sec_send_double[XDOWN];
		break;
	case sec_send_yup:
		result = g_bgq_sec_send_double[YUP];
		break;
	case sec_send_ydown:
		result = g_bgq_sec_send_double[YDOWN];
		break;
	case sec_send_zup:
		result = g_bgq_sec_send_double[ZUP];
		break;
	case sec_send_zdown:
		result = g_bgq_sec_send_double[ZDOWN];
		break;
	case sec_recv_tup:
		result = g_bgq_sec_recv_double[TUP];
		break;
	case sec_recv_tdown:
		result = g_bgq_sec_recv_double[TDOWN];
		break;
	case sec_recv_xup:
		result = g_bgq_sec_recv_double[XUP];
		break;
	case sec_recv_xdown:
		result = g_bgq_sec_recv_double[XDOWN];
		break;
	case sec_recv_yup:
		result = g_bgq_sec_recv_double[YUP];
		break;
	case sec_recv_ydown:
		result = g_bgq_sec_recv_double[YDOWN];
		break;
	case sec_recv_zup:
		result = g_bgq_sec_recv_double[ZUP];
		break;
	case sec_recv_zdown:
		result = g_bgq_sec_recv_double[ZDOWN];
		break;
	case sec_temp_tup:
		result = g_bgq_sec_temp_tup_double;
		break;
	case sec_temp_tdown:
		result = g_bgq_sec_temp_tdown_double;
		break;
	case sec_vrecv_tdown:
	case sec_vrecv_tup:
		return NULL;
	default:
		assert(!"Section does not physically exist");
		UNREACHABLE
	}

	return result;
}


bgq_weyl_vec_float *bgq_section_baseptr_float(bgq_weylfield_controlblock *field, bgq_weylfield_section section) {
	assert(field);
	//assert(field->isInitialized);
	bgq_weyl_vec_float *result = NULL;

	switch (section) {
	case sec_surface:
		if (PHYSICAL_SURFACE==0)
			return NULL;
		result = (bgq_weyl_vec_float*)&field->sec_collapsed_float[bgq_surface2collapsed(0)];
		break;
	case sec_body:
		if (PHYSICAL_BODY==0)
			return NULL;
		result = (bgq_weyl_vec_float*)&field->sec_collapsed_float[bgq_body2collapsed(0)];
		break;
	case sec_send_tup:
		result = g_bgq_sec_send_float[TUP];
		break;
	case sec_send_tdown:
		result = g_bgq_sec_send_float[TDOWN];
		break;
	case sec_send_xup:
		result = g_bgq_sec_send_float[XUP];
		break;
	case sec_send_xdown:
		result = g_bgq_sec_send_float[XDOWN];
		break;
	case sec_send_yup:
		result = g_bgq_sec_send_float[YUP];
		break;
	case sec_send_ydown:
		result = g_bgq_sec_send_float[YDOWN];
		break;
	case sec_send_zup:
		result = g_bgq_sec_send_float[ZUP];
		break;
	case sec_send_zdown:
		result = g_bgq_sec_send_float[ZDOWN];
		break;
	case sec_recv_tup:
		result = g_bgq_sec_recv_float[TUP];
		break;
	case sec_recv_tdown:
		result = g_bgq_sec_recv_float[TDOWN];
		break;
	case sec_recv_xup:
		result = g_bgq_sec_recv_float[XUP];
		break;
	case sec_recv_xdown:
		result = g_bgq_sec_recv_float[XDOWN];
		break;
	case sec_recv_yup:
		result = g_bgq_sec_recv_float[YUP];
		break;
	case sec_recv_ydown:
		result = g_bgq_sec_recv_float[YDOWN];
		break;
	case sec_recv_zup:
		result = g_bgq_sec_recv_float[ZUP];
		break;
	case sec_recv_zdown:
		result = g_bgq_sec_recv_float[ZDOWN];
		break;
	case sec_temp_tup:
		result = g_bgq_sec_temp_tup_float;
		break;
	case sec_temp_tdown:
		result = g_bgq_sec_temp_tdown_float;
		break;
	case sec_vrecv_tdown:
	case sec_vrecv_tup:
		return NULL;
	default:
		assert(!"Section does not physically exist");
		UNREACHABLE
	}

	return result;
}


static bgq_spinorfield_layout bgq_spinorfield_bestLayout(bgq_weylfield_controlblock *field) {
	assert(field);
	//assert(field->isInitialized);

	if (field->has_fulllayout_double) {
		return ly_full_double;
	} else if (field->has_weyllayout_double) {
		return ly_weyl_double;
	} else if (field->has_fulllayout_float) {
		return ly_full_float;
	} else if (field->has_weyllayout_float) {
		return ly_weyl_float;
	} else if (field->has_legacy) {
		return ly_legacy;
	}
	return ly_none;
}


static inline void bgq_copyToLegacy_worker(void *arg_untyped, size_t tid, size_t threads, bool weyllayout, bool sloppy, bool mul, bool isLegacy) {
	bgq_copyToLegacy_workload *arg = arg_untyped;
	bool isOdd = arg->isOdd;
	bgq_weylfield_controlblock *field = arg->field;
	spinor *legacyField = arg->target;

	assert(field->isOdd == isOdd);

	const size_t workload = PHYSICAL_VOLUME;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);

	bgq_spinorfield_streamSpinor(field, isOdd, begin, weyllayout, sloppy, mul, isLegacy);
	for (ucoord ic = begin; ic<end; ic+=1) {
#ifndef NDEBUG
		ucoord t1 = bgq_collapsed2t1(isOdd, ic);
		ucoord t2 = bgq_collapsed2t2(isOdd, ic);
		ucoord x = bgq_collapsed2x(isOdd, ic);
		ucoord y = bgq_collapsed2y(isOdd, ic);
		ucoord z = bgq_collapsed2z(isOdd, ic);
#endif

		bgq_su3_spinor_decl(spinor);
		bgq_spinorfield_readSpinor(&spinor, field, isOdd, ic, weyllayout, sloppy, mul, isLegacy);
				bgq_spinorqpx_expect(spinor,t1,t2,x,y,z);
		bgq_spinorfield_prefetchNextSpinor(field, isOdd, ic, weyllayout, sloppy, mul, isLegacy);

		int eosub1 = bgq_collapsed2eosub(isOdd, ic, 0);
		assert(bgq_eosub2collapsed(isOdd, eosub1) == ic);
		spinor *addr1 = &legacyField[eosub1];
		bgq_st2a_double(spinor_v0_c0, 0, addr1);
		bgq_qvstfcduxa(spinor_v0_c1, addr1, 16);
		bgq_qvstfcduxa(spinor_v0_c2, addr1, 16);
		bgq_qvstfcduxa(spinor_v1_c0, addr1, 16);
		bgq_qvstfcduxa(spinor_v1_c1, addr1, 16);
		bgq_qvstfcduxa(spinor_v1_c2, addr1, 16);
		bgq_qvstfcduxa(spinor_v2_c0, addr1, 16);
		bgq_qvstfcduxa(spinor_v2_c1, addr1, 16);
		bgq_qvstfcduxa(spinor_v2_c2, addr1, 16);
		bgq_qvstfcduxa(spinor_v3_c0, addr1, 16);
		bgq_qvstfcduxa(spinor_v3_c1, addr1, 16);
		bgq_qvstfcduxa(spinor_v3_c2, addr1, 16);
				bgq_spinorlegacy_expect(&legacyField[eosub1], t1, x, y, z);

		int eosub2 = bgq_collapsed2eosub(isOdd, ic, 1);
		assert(bgq_eosub2collapsed(isOdd, eosub2) == ic);
		spinor *addr2 = &legacyField[eosub2];
		bgq_vector4double_decl(right0);
		bgq_vector4double_decl(right1);
		bgq_vector4double_decl(right2);
		bgq_vector4double_decl(right3);
		bgq_vector4double_decl(right4);
		bgq_vector4double_decl(right5);
		bgq_rmerge(right0, spinor_v0_c0, spinor_v0_c1);
		bgq_rmerge(right1, spinor_v0_c2, spinor_v1_c0);
		bgq_rmerge(right2, spinor_v1_c1, spinor_v1_c2);
		bgq_rmerge(right3, spinor_v2_c0, spinor_v2_c1);
		bgq_rmerge(right4, spinor_v2_c2, spinor_v3_c0);
		bgq_rmerge(right5, spinor_v3_c1, spinor_v3_c2);
		bgq_sta_double(right0, 0, addr2);
		bgq_qvstfduxa(right1, addr2, 32);
		bgq_qvstfduxa(right2, addr2, 32);
		bgq_qvstfduxa(right3, addr2, 32);
		bgq_qvstfduxa(right4, addr2, 32);
		bgq_qvstfduxa(right5, addr2, 32);
				bgq_spinorlegacy_expect(&legacyField[eosub2], t2, x, y, z);
	}
}

BGQ_SPINORFIELD_GENWORKER(bgq_copyToLegacy_worker)


void bgq_spinorfield_prepareWrite(bgq_weylfield_controlblock *field, tristate isOdd, bgq_spinorfield_layout layout, bool preserveData) {
	assert(field);

	if (preserveData) {
		if (isOdd==tri_unknown) {
			isOdd = field->isOdd;
		} else if (field->isOdd==tri_unknown) {
		} else {
			assert(isOdd==field->isOdd);
		}
	}

	if (layout & ly_weyl) {
		// Ensure any communication has finished before messing up with fields
		// Empty communication buffer so they can be reused
		bgq_comm_wait();
	}
	else if (field->pendingDatamove) {
		// This is a bit strange; actually it means that the result of a HoppingMatrix call has not been used
		bgq_comm_wait();
	}

#ifndef NDEBUG
	// Some debugging code in workers (especially the datamove inside bgq_comm_wait() called above) needs to know whether the float or double field is currently active
	// But exactly this we are going to change here, so in debug mode, wait for workers here
	// Pay attention because this might mean that there could be race condition that exist in the release version, but not in debug version
	// Alternative: In debug mode, use different malloc for float and double fields
	bgq_master_sync();
#endif

	bgq_spinorfield_enableLayout(field, isOdd, layout, true, preserveData);
}


void bgq_spinorfield_prepareReadWrite(bgq_weylfield_controlblock *field, tristate isOdd, bgq_spinorfield_layout layout) {
	isOdd = tristate_combine(isOdd, field->isOdd);
	bgq_spinorfield_prepareRead(field, isOdd,  layout&ly_weyl, (layout!=ly_legacy) && !(layout&ly_sloppy), layout&ly_sloppy, layout&ly_mul, layout==ly_legacy);
	bgq_spinorfield_prepareWrite(field, isOdd, layout, true);
}


bgq_spinorfield_layout bgq_spinorfield_prepareRead(bgq_weylfield_controlblock *field, tristate isOdd, bool acceptWeyl, bool acceptDouble, bool acceptFloat, bool acceptMul, bool acceptLegacy) {
	assert(field);
	assert(field->has_fulllayout_double || field->has_fulllayout_float || field->has_weyllayout_double || field->has_weyllayout_float || field->has_legacy); // There must be some data to read
	assert(acceptDouble || acceptFloat || acceptLegacy); // Accept at least something

	if (isOdd == field->isOdd) {
		// Matching oddness, or unknown
	} else if (isOdd==tri_unknown) {
		isOdd = field->isOdd;
	} else if (field->isOdd==tri_unknown) {
		field->isOdd = isOdd;
	} else {
		// Mismatching oddness
		master_error(1, "Oddness mismatch");
	}

	bool actionRewrite = false;
	bgq_spinorfield_layout layout;
	if (field->has_fulllayout_double && acceptDouble) {
		layout = ly_full_double;
	} else if (field->has_weyllayout_double && acceptWeyl && acceptDouble) {
		layout = ly_weyl_double;
	} else if (field->has_fulllayout_float  && acceptFloat) {
		layout = ly_full_float;
	} else if (field->has_weyllayout_float && acceptWeyl && acceptFloat) {
		layout = ly_weyl_float;
	} else if (field->has_legacy && acceptLegacy) {
		layout = ly_legacy;
	} else {
		actionRewrite = true;
		layout = bgq_spinorfield_bestLayout(field);
	}
	assert(layout!=ly_none);

	if (field->pendingDatamove) {
		// 4, Wait for data to be received
		bgq_comm_wait();
	}

	bgq_spinorfield_layout result = -1;
	if (actionRewrite) {
		if (isOdd==tri_unknown) {
			master_error(1, "ERROR: For rewriting, we really need to know the oddness of the field from somewhere\n");
		}
		if (acceptDouble) {
			result = ly_full_double;
			bgq_worker_func worker = g_bgq_spinorfield_rewrite_worker_double_list[layout];
			assert(worker);

			if ((layout==ly_full_float) && ((void*)field->sec_fullspinor_double==(void*)field->sec_fullspinor_float)) {
				// Without intervention, we're going to overwrite the data while we are reading it
				// Solution: assign a new memory area
				// Note that field->sec_fullspinor_float has allocated twice as much memory, so we are wasting some space here
				field->sec_fullspinor_double = malloc_aligned(LOCAL_VOLUME/PHYSICAL_LP * sizeof(*field->sec_fullspinor_double), BGQ_ALIGNMENT_L2);
#ifndef NVALGRIND
				VALGRIND_CREATE_MEMPOOL(field->sec_fullspinor_double, 0, false);
				VALGRIND_MEMPOOL_ALLOC(field->sec_fullspinor_double, field->sec_fullspinor_double, PHYSICAL_VOLUME * sizeof(*field->sec_fullspinor_double));
#endif
			}
			bgq_spinorfield_enableLayout(field, isOdd, result, false, false);

			bgq_master_sync();
			static bgq_spinorfield_rewrite_work work;
			work.isOdd = isOdd;
			work.field = field;
			bgq_master_call(worker, &work);
		} else if (acceptFloat) {
			result = ly_full_float;
			bgq_worker_func worker = g_bgq_spinorfield_rewrite_worker_float_list[layout];
			assert(worker);

			if ((layout==ly_full_double) && ((void*)field->sec_fullspinor_double==(void*)field->sec_fullspinor_float)) {
				// Without intervention, we're going to overwrite the data while we are reading it
				// Solution: assign a new memory area
				field->sec_fullspinor_float = malloc_aligned(LOCAL_VOLUME/PHYSICAL_LP * sizeof(*field->sec_fullspinor_float), BGQ_ALIGNMENT_L2);
#ifndef NVALGRIND
				VALGRIND_CREATE_MEMPOOL(field->sec_fullspinor_float, 0, false);
				VALGRIND_MEMPOOL_ALLOC(field->sec_fullspinor_float, field->sec_fullspinor_float, PHYSICAL_VOLUME * sizeof(*field->sec_fullspinor_float));
#endif
			}
			bgq_spinorfield_enableLayout(field, isOdd, result, false, false);

			bgq_master_sync();
			static bgq_spinorfield_rewrite_work work;
			work.isOdd = isOdd;
			work.field = field;
			bgq_master_call(worker, &work);
		} else if (acceptLegacy) {
			result = ly_legacy;
			// We have to copy existing data
			assert(layout != ly_legacy);
			//master_print("Translation to ly_legacy\n");

			if (field->pendingDatamove) {
				bgq_comm_wait();
			}

			//bgq_spinorfield_layout layout = bgq_spinorfield_bestLayout(field);
			bgq_spinorfield_enableLayout(field, isOdd, ly_legacy, false, false);
			bgq_worker_func worker = g_bgq_copyToLegacy_worker_list[layout];
			assert(worker);

			bgq_master_sync();
			static bgq_copyToLegacy_workload work;
			work.isOdd = isOdd;
			work.field = field;
			work.target = field->legacy_field;
			bgq_master_call(worker, &work);

#ifndef NDEBUG
			 double diff = bgq_spinorfield_legacy_compare(isOdd,field, layout, field->legacy_field, false);
			 assert(diff==0);
#endif
		} else {
			assert(!"You ain't accept anything we can convert to!");
		}
	} else {
		result = layout;
	}

	assert(result!=ly_none);
	assert(acceptWeyl || !(result&ly_weyl));
	assert(acceptDouble || (result&ly_sloppy) || (acceptLegacy && result==ly_legacy));
	assert(acceptFloat || !(result&ly_sloppy) || (acceptLegacy && result==ly_legacy));
	assert(acceptMul || !(result&ly_mul));
	assert(acceptLegacy || (result!=ly_legacy));

	return result;
}


size_t bgq_pointer2offset_raw(bgq_weylfield_controlblock *field, void *ptr, bool check) {
	if (!field) {
		if (check) {
			assert(!"No field passed");
		}
		return -1;
	}

	for (bgq_weylfield_section sec = 0; sec < sec_end; sec+=1) {
		size_t secsize_double = bgq_section_size(sec);
		size_t secsize_float = secsize_double/2;

		if (field->has_weyllayout_float && field->has_weyllayout_double && ((void*)field->sec_collapsed_double == (void*)field->sec_collapsed_float)) {
			assert(!"Double and float pointers indistinguishable");
		}

		if (field->has_weyllayout_float) {
			bgq_weyl_vec_float *baseptr = bgq_section_baseptr_float(field, sec);
			if ((uint8_t*)baseptr <= (uint8_t*)ptr && (uint8_t*)ptr < (uint8_t*)baseptr+secsize_float) {
				size_t baseoffset = bgq_weyl_section_offset(sec);
				size_t result = baseoffset + 2*((uint8_t*)ptr - (uint8_t*)baseptr);
				assert(result % sizeof(bgq_weyl_vec_float) == 0);
				assert(result < bgq_weyl_section_offset(sec_end));
				return result;
			}
		}

		if (field->has_weyllayout_double) {
			bgq_weyl_vec_double *baseptr = bgq_section_baseptr_double(field, sec);
			if ((uint8_t*)baseptr <= (uint8_t*)ptr && (uint8_t*)ptr < (uint8_t*)baseptr+secsize_double) {
				size_t baseoffset = bgq_weyl_section_offset(sec);
				size_t result = baseoffset + ((uint8_t*)ptr - (uint8_t*)baseptr);
				assert(result % sizeof(bgq_weyl_vec_double) == 0);
				assert(result < bgq_weyl_section_offset(sec_end));
				return result;
			}
		}
	}
	if (check) {
		assert(!"Pointer to non-field location");
	}
	return -1;
}


size_t bgq_pointer2offset(bgq_weylfield_controlblock *field, void *ptr) {
	return bgq_pointer2offset_raw(field, ptr, true);
}





static bool g_bgq_spinorfields_initialized = false;


static bgq_weylfield_collection *g_bgq_spinorfield_collection_first = NULL;
static bgq_weylfield_collection *g_bgq_spinorfield_collection_last = NULL;

static bgq_weylfield_collection *g_bgq_spinorfield_collection_unused = NULL;


bgq_weylfield_collection *bgq_spinorfields_allocate(size_t count, spinor *legacyFields, size_t fieldLength) {
	bgq_spinorfields_init();
	assert(fieldLength >= VOLUME/2);

	bgq_weylfield_collection *result;
	size_t nInitialized = 0;
	if (g_bgq_spinorfield_collection_unused) {
		result = g_bgq_spinorfield_collection_unused;
		g_bgq_spinorfield_collection_unused = result->next;
		if (g_bgq_spinorfield_collection_unused) {
			g_bgq_spinorfield_collection_unused->prev = NULL;
		}

		nInitialized = result->count;
		if (result->count < count) {
			result = realloc(result, sizeof(*result) + count * sizeof(*result->controlblocks));
		} else {
			count = result->count;
		}
	} else {
		result = malloc(sizeof(*result) + count * sizeof(*result->controlblocks));
	}

	result->count = count;
	result->legacy_base = legacyFields;
	result->prev = NULL;
	result->next = NULL;

	size_t fieldsize = /*VOLUMEPLUSRAND/2*/ fieldLength * sizeof(spinor);
	result->fieldsize = fieldsize;

#ifndef NVALGRIND
	VALGRIND_CREATE_MEMPOOL(legacyFields, 0, false);
#endif

	for (size_t i = 0; i < count; i += 1) {
		bgq_weylfield_controlblock *field = &result->controlblocks[i];

		field->has_legacy = true;
		field->has_fulllayout_double = false;
		field->has_fulllayout_float = false;
		field->has_weyllayout_double = false;
		field->has_weyllayout_float = false;

		if (i < nInitialized) {
		} else {
			field->sec_collapsed_double = NULL;
			field->sec_collapsed_float = NULL;
			field->sec_fullspinor_double = NULL;
			field->sec_fullspinor_float = NULL;

			for (size_t isOdd = 0; isOdd < PHYSICAL_LP; isOdd+=1) {
				field->sendptr_double[isOdd] = NULL;
				field->sendptr_float[isOdd] = NULL;
				for (size_t d = 0; d < PHYSICAL_LD; d+=1) {
					field->consptr_double[isOdd][d] = NULL;
					field->consptr_float[isOdd][d] = NULL;
				}
			}
		}

		field->collectionBase = result;
		field->isOdd = tri_unknown;

		if (legacyFields) {
			field->legacy_field = (spinor*)((uint8_t*)legacyFields + i*fieldsize);
			assert((uintptr_t)field->legacy_field % 32 == 0);
#ifndef NVALGRIND
			VALGRIND_MEMPOOL_ALLOC(legacyFields, field->legacy_field, VOLUMEPLUSRAND/2*sizeof(*field->legacy_field));
#endif
#ifndef NDEBUG
			memset(field->legacy_field, 0xFF, fieldsize);
#endif
		} else {
			field->legacy_field = NULL;
		}

		field->pendingDatamove = false;
	}

	if (g_bgq_spinorfield_collection_last) {
		assert(g_bgq_spinorfield_collection_last->next==NULL);
		g_bgq_spinorfield_collection_last->next = result;
		result->prev = g_bgq_spinorfield_collection_last;
		g_bgq_spinorfield_collection_last = result;
	} else {
		g_bgq_spinorfield_collection_first = result;
		g_bgq_spinorfield_collection_last = result;
	}
	assert(!g_bgq_spinorfield_collection_first || !g_bgq_spinorfield_collection_first->prev);

	return result;
}


void bgq_spinorfields_free(bgq_weylfield_collection *collection) {
	assert(collection);

	bgq_weylfield_collection *prev = collection->prev;
	bgq_weylfield_collection *next = collection->next;

	if (prev) {
		prev->next = (prev==next) ? NULL : next;
	}
	if (next) {
		next->prev = (prev==next) ? NULL : prev;
	}

	if (collection==g_bgq_spinorfield_collection_first) {
		assert(!prev);
		g_bgq_spinorfield_collection_first = next;
	}

	if (collection==g_bgq_spinorfield_collection_last) {
		assert(!next);
		g_bgq_spinorfield_collection_last = prev;
	}

	collection->prev = NULL;
	collection->next = NULL;

	if (g_bgq_spinorfield_collection_unused) {
		g_bgq_spinorfield_collection_unused->prev = collection;
		collection->next = g_bgq_spinorfield_collection_unused;
		g_bgq_spinorfield_collection_unused = collection;
	} else {
		g_bgq_spinorfield_collection_unused = collection;
	}


#ifndef NVALGRIND
	VALGRIND_DESTROY_MEMPOOL(collection->legacy_base);
#endif
}


void bgq_spinorfields_init(void) {
	bgq_indices_init();
	if (g_bgq_spinorfields_initialized)
		return;
	g_bgq_spinorfields_initialized = true;

	//bgq_weylfield_collection *collection = bgq_spinorfields_allocate(std_count, g_spinor_field[0], VOLUMEPLUSRAND/2);
	//g_bgq_spinorfields = &collection->controlblocks[0];
}


bgq_weylfield_controlblock *bgq_translate_spinorfield(const spinor *legacy_field) {
	uint8_t *legacyField = (uint8_t*)legacy_field;

	bgq_weylfield_collection *collection = g_bgq_spinorfield_collection_first;
	while (collection) {
		size_t count = collection->count;
		size_t fieldsize = collection->fieldsize;// VOLUMEPLUSRAND/2 * sizeof(spinor);
		uint8_t *legacyBase = (uint8_t*)collection->legacy_base;

		if (legacyField < legacyBase) {
			collection = collection->next;
			continue;
		}
		ptrdiff_t offset = legacyField - legacyBase;
		if (offset >= count * fieldsize) {
			printf("ciao %p %p %d %d\n",legacyField,legacyBase,count,fieldsize);
			collection = collection->next;
			continue;
		}

		assert(offset % fieldsize == 0);
		size_t index = offset / fieldsize;

		bgq_weylfield_controlblock *result = &collection->controlblocks[index];
		assert(result->legacy_field == (spinor*)legacyField);
		assert(result!=NULL);
		return result;
	}

	assert(false);
	return NULL;
}


void spinorfield_enable(const spinor *legacyField, int read, int write) {
	bgq_weylfield_controlblock *field = bgq_translate_spinorfield(legacyField);
	assert(read || write);

	if (read) {
		bgq_spinorfield_prepareRead(field, tri_unknown, false, false, false, false, true);
	}

	if (write) {
		bgq_spinorfield_prepareWrite(field, tri_unknown, ly_legacy, read);
	}
}


void spinorfield_propagateOddness(const spinor *targetLegacyField, const spinor *sourceLegacyField) {
	assert(targetLegacyField!=sourceLegacyField); // Probably field equality has not been chacked when calling spinorfield_enable
	if (targetLegacyField==sourceLegacyField)
		return;

	bgq_weylfield_controlblock *targetField = bgq_translate_spinorfield(targetLegacyField);
	bgq_weylfield_controlblock *sourceField = bgq_translate_spinorfield(sourceLegacyField);

	if (targetField->isOdd==tri_unknown) {
		targetField->isOdd = sourceField->isOdd;
	} else if (sourceField->isOdd==tri_unknown) {
	} else {
		assert(sourceField->isOdd==targetField->isOdd);
	}
}


void spinorfield_propagateInvertedOddness(const spinor *targetLegacyField, const spinor *sourceLegacyField) {
	if (targetLegacyField==sourceLegacyField)
		return;

	bgq_weylfield_controlblock *targetField = bgq_translate_spinorfield(targetLegacyField);
	bgq_weylfield_controlblock *sourceField = bgq_translate_spinorfield(sourceLegacyField);

	if (sourceField->isOdd==tri_unknown) {
	} else if (targetField->isOdd==tri_unknown) {
		targetField->isOdd = !sourceField->isOdd;
	} else {
		assert(sourceField->isOdd==!targetField->isOdd);
	}
}


void bgq_legacy_markcoords_raw(bool isOdd, spinor *legacyField) {
	spinorfield_enable(legacyField, false, true);

	for (int eosub = 0; eosub < VOLUME/2; eosub += 1) {
		int ioff = isOdd ? (VOLUME+RAND)/2 : 0;
		int eo = eosub + ioff;
		int lexic = g_eo2lexic[eo];
		int t = g_coord[lexic][0];
		int x = g_coord[lexic][1];
		int y = g_coord[lexic][2];
		int z = g_coord[lexic][3];

		legacyField[eosub].s0.c0 = t;
		legacyField[eosub].s0.c1 = x;
		legacyField[eosub].s0.c2 = y;
		legacyField[eosub].s1.c0 = z;
		legacyField[eosub].s1.c1 = 0;
		legacyField[eosub].s1.c2 = 0;
		legacyField[eosub].s2.c0 = 0;
		legacyField[eosub].s2.c1 = 0;
		legacyField[eosub].s2.c2 = 0;
		legacyField[eosub].s3.c0 = 0;
		legacyField[eosub].s3.c1 = 0;
		legacyField[eosub].s3.c2 = 0;
	}
}


void bgq_spinorfield_zero(bgq_weylfield_controlblock *field, tristate isOdd) {
	assert(field);

	bgq_spinorfield_enableLayout(field, isOdd, ly_legacy, false, false); // Just to set oddness
	bgq_master_memzero(field->legacy_field, VOLUME/2 * sizeof(*field->legacy_field));
	field->has_legacy = true;

	//bgq_spinorfield_enableLayout(field, isOdd, ly_full_double, false, false);
	if (field->sec_fullspinor_double) {
		bgq_master_memzero(field->sec_fullspinor_double, LOCAL_VOLUME/PHYSICAL_LP * sizeof(*field->sec_fullspinor_double));
		field->has_fulllayout_double = true;
	}

	if (field->sec_fullspinor_float) {
		if ((void*)field->sec_fullspinor_float!=(void*)field->sec_fullspinor_double)
			bgq_master_memzero(field->sec_fullspinor_float, LOCAL_VOLUME/PHYSICAL_LP * sizeof(*field->sec_fullspinor_float));
		field->has_fulllayout_float = true;
	}

	if (field->sec_collapsed_double) {
		bgq_master_memzero(field->sec_collapsed_double, PHYSICAL_VOLUME * sizeof(*field->sec_collapsed_double));
		field->has_weyllayout_double = true;
	}

	if (field->sec_collapsed_float) {
		if ((void*)field->sec_collapsed_float!=(void*)field->sec_collapsed_double)
			bgq_master_memzero(field->sec_collapsed_float, PHYSICAL_VOLUME * sizeof(*field->sec_collapsed_float));
		field->has_weyllayout_float = true;
	}
}


void bgq_spinorfield_annotateOddness(bgq_weylfield_controlblock *field, bool isOdd) {
	assert(field);

	if (field->isOdd == tri_unknown) {
		field->isOdd = isOdd;
	} else {
		if (field->isOdd != isOdd) {
			master_error(1, "Mismatch in oddness");
		}
	}
}


void spinorfield_setOddness(const spinor *field, int isOdd) {
	bgq_weylfield_controlblock *sfield = bgq_translate_spinorfield(field);
	bgq_spinorfield_annotateOddness(sfield, isOdd);
}


void spinorfield_linalg_u(const spinor *legacyField_inout) {
	bgq_weylfield_controlblock *field_inout = bgq_translate_spinorfield(legacyField_inout);

	bgq_spinorfield_prepareReadWrite(field_inout, tri_unknown, ly_legacy);
}

/* write read read */
void spinorfield_linalg_wrr(const spinor *legacyField_out, const spinor *legacyField_in1, const spinor *legacyField_in2) {
	bgq_weylfield_controlblock *field_out = bgq_translate_spinorfield(legacyField_out);
	bgq_weylfield_controlblock *field_in1 = bgq_translate_spinorfield(legacyField_in1);
	bgq_weylfield_controlblock *field_in2 = bgq_translate_spinorfield(legacyField_in2);

	tristate isOdd = tristate_combine(field_in1->isOdd, field_in2->isOdd);
	bgq_spinorfield_prepareRead(field_in1, isOdd, false, false, false, false, true);
	bgq_spinorfield_prepareRead(field_in2, isOdd, false, false, false, false, true);
	bgq_spinorfield_prepareWrite(field_out, isOdd, ly_legacy, (field_out==field_in1) || (field_out==field_in2));
}

/* reUse read read*/
void spinorfield_linalg_urr(const spinor *legacyField_inout, const spinor *legacyField_in1, const spinor *legacyField_in2) {
	bgq_weylfield_controlblock *field_inout = bgq_translate_spinorfield(legacyField_inout);
	bgq_weylfield_controlblock *field_in1 = bgq_translate_spinorfield(legacyField_in1);
	bgq_weylfield_controlblock *field_in2 = bgq_translate_spinorfield(legacyField_in2);

	tristate isOdd = tristate_combine(field_inout->isOdd, tristate_combine(field_in1->isOdd, field_in2->isOdd));
	bgq_spinorfield_prepareRead(field_in1, isOdd, false, false, false, false, true);
	bgq_spinorfield_prepareRead(field_in2, isOdd, false, false, false, false, true);
	bgq_spinorfield_prepareReadWrite(field_inout, isOdd, ly_legacy);
}

/* reUse read */
void spinorfield_linalg_ur(const spinor *legacyField_inout, const spinor *legacyField_in) {
	bgq_weylfield_controlblock *field_inout = bgq_translate_spinorfield(legacyField_inout);
	bgq_weylfield_controlblock *field_in = bgq_translate_spinorfield(legacyField_in);

	tristate isOdd = tristate_combine(field_inout->isOdd, field_in->isOdd);
	bgq_spinorfield_prepareRead(field_in, isOdd, false, false, false, false, true);
	bgq_spinorfield_prepareReadWrite(field_inout, isOdd, ly_legacy);
}

/* reUse read read read */
void spinorfield_linalg_urrr(const spinor *legacyField_inout, const spinor *legacyField_in1, const spinor *legacyField_in2, const spinor *legacyField_in3) {
	bgq_weylfield_controlblock *field_inout = bgq_translate_spinorfield(legacyField_inout);
	bgq_weylfield_controlblock *field_in1 = bgq_translate_spinorfield(legacyField_in1);
	bgq_weylfield_controlblock *field_in2 = bgq_translate_spinorfield(legacyField_in2);
	bgq_weylfield_controlblock *field_in3 = bgq_translate_spinorfield(legacyField_in3);

	tristate isOdd = tristate_combine(field_inout->isOdd, tristate_combine(tristate_combine(field_in1->isOdd, field_in2->isOdd), field_in3->isOdd));
	bgq_spinorfield_prepareRead(field_in1, isOdd, false, false, false, false, true);
	bgq_spinorfield_prepareRead(field_in2, isOdd, false, false, false, false, true);
	bgq_spinorfield_prepareRead(field_in3, isOdd, false, false, false, false, true);
	bgq_spinorfield_prepareReadWrite(field_inout, isOdd, ly_legacy);
}

/* write read */
void spinorfield_linalg_wr(const spinor *legacyField_out, const spinor *legacyField_in) {
	bgq_weylfield_controlblock *field_out = bgq_translate_spinorfield(legacyField_out);
	bgq_weylfield_controlblock *field_in = bgq_translate_spinorfield(legacyField_in);

	tristate isOdd = field_in->isOdd;
	bgq_spinorfield_prepareRead(field_in, isOdd, false, false, false, false, true);
	bgq_spinorfield_prepareWrite(field_out, isOdd, ly_legacy, field_out==field_in);
}

/* write read */
void spinorfield_linalg_rr(const spinor *field_in1, const spinor *field_in2) {
	spinorfield_enable(field_in1, 1, 0);
	spinorfield_enable(field_in2, 1, 0);
}

/* write write */
void spinorfield_linalg_ww(const spinor *field_out1, const spinor *field_out2) {
	spinorfield_enable(field_out1, 0, 1);
	spinorfield_enable(field_out2, 0, 1);
}

/* read */
void spinorfield_linalg_r(const spinor *field_in) {
	spinorfield_enable(field_in, 1, 0);
}

/* write read */
void spinorfield_linalg_wr_invert(const spinor *legacyField_out, const spinor *legacyField_in) {
	assert(legacyField_out != legacyField_in);
	bgq_weylfield_controlblock *field_out = bgq_translate_spinorfield(legacyField_out);
	bgq_weylfield_controlblock *field_in = bgq_translate_spinorfield(legacyField_in);

	tristate isOdd = field_in->isOdd;
	bgq_spinorfield_prepareRead(field_in, isOdd, false, false, false, false, true);
	bgq_spinorfield_prepareWrite(field_out, tristate_invert(isOdd), ly_legacy, false);
}

/* write read read */
void spinorfield_linalg_wrr_invert(const spinor *legacyField_out, const spinor *legacyField_in1, const spinor *legacyField_in2) {
	assert(legacyField_out != legacyField_in1);
	assert(legacyField_out != legacyField_in2);
	bgq_weylfield_controlblock *field_out = bgq_translate_spinorfield(legacyField_out);
	bgq_weylfield_controlblock *field_in1 = bgq_translate_spinorfield(legacyField_in1);
	bgq_weylfield_controlblock *field_in2 = bgq_translate_spinorfield(legacyField_in2);

	tristate isOdd = tristate_combine(tristate_invert(field_in1->isOdd), field_in2->isOdd);
	bgq_spinorfield_prepareRead(field_in1, tristate_invert(isOdd), false, false, false, false, true);
	bgq_spinorfield_prepareRead(field_in2, isOdd, false, false, false, false, true);
	bgq_spinorfield_prepareWrite(field_out, tristate_invert(isOdd), ly_legacy, false);
}


