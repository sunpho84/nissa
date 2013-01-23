/*
 * bgq_field.c
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */


#define BGQ_FIELD_C_
#include "bgq_field.h"

#include "bgq_HoppingMatrix.h"
#include "bgq_dispatch.h"
#include "bgq_qpx.h"
#include "bgq_comm.h"




static bgq_weylfield_section bgq_HoppingMatrix_init_source_sectionof_local(size_t t, size_t x, size_t y, size_t z, bgq_direction d) {
	// Check whether moving out of the local lattice
	if (HALO_T && (t==LOCAL_LT-1) && (d==TUP)) {
		if (BGQ_UNVECTORIZE || !COMM_T)
			return sec_temp_tup;
		else
			return sec_send_tup;
	}
	if (HALO_T && (t==0) && (d==TDOWN)) {
		if (BGQ_UNVECTORIZE || !COMM_T)
			return sec_temp_tdown;
		else
			return sec_send_tdown;
	}
	if (HALO_X && (x==LOCAL_LX-1) && (d==XUP)) {
		return sec_send_xup;
	}
	if (HALO_X && (x==0) && (d==XDOWN)) {
		return sec_send_xdown;
	}
	if (HALO_Y && (y==LOCAL_LY-1) && (d==YUP)) {
		return sec_send_yup;
	}
	if (HALO_Y && (y==0) && (d==YDOWN)) {
		return sec_send_ydown;
	}
	if (HALO_Z && (z==LOCAL_LZ-1) && (d==ZUP)) {
		return sec_send_zup;
	}
	if (HALO_Z && (z==0) && (d==ZDOWN)) {
		return sec_send_zdown;
	}

	// Stay within the local lattice, determine whether the target is on the surface
	bgq_direction_move_local(&t, &x, &y, &z, d);
	bool isSurface = bgq_local2isSurface(t, x, y, z);
	return isSurface ? sec_surface : sec_body;
}


static bgq_weylfield_section bgq_HoppingMatrix_init_source_sectionof_physical(bool isOdd, size_t tv, size_t x, size_t y, size_t z, bgq_direction d) {
	size_t t1 = bgq_physical2t1(isOdd, tv, x, y, z);
	bgq_weylfield_section sec1 = bgq_HoppingMatrix_init_source_sectionof_local(t1,x,y,z,d);
	size_t t2 = bgq_physical2t2(isOdd, tv, x, y, z);
	bgq_weylfield_section sec2 = bgq_HoppingMatrix_init_source_sectionof_local(t2,x,y,z,d);

	if (sec1 == sec2) {
		return sec1;
	} else if (((sec1 == sec_surface) || (sec1 == sec_body)) && ((sec2 == sec_send_tup) || (sec2 == sec_temp_tup))) {
		// Happens because t1 is not strictly at the border
		// We give t2 priority such that we save both into the buffer
		return sec2;
	} else if (((sec2 == sec_surface) || (sec2 == sec_body)) && ((sec1 == sec_send_tdown) || (sec1 == sec_temp_tdown))) {
		return sec1;
	} else {
		assert(!"Unknown case");
		UNREACHABLE
		return -1;
	}
}


static size_t bgq_offset_send2recv(size_t offset_send) {
	bgq_weylfield_section sec_send = bgq_sectionOfOffset(offset_send);
	if ((sec_send==sec_surface)||(sec_send==sec_body))
		return offset_send;
	if (!COMM_T && ((sec_send==sec_recv_tdown) || (sec_send==sec_send_tup)))
		return offset_send;
	if (!COMM_T && ((sec_send==sec_recv_tup) || (sec_send==sec_send_tdown)))
		return offset_send;
	bgq_weylfield_section sec_recv = bgq_section_commbuftran(sec_send, false);
	assert(sec_recv!=sec_send);

	size_t offset_recv = offset_send - bgq_weyl_section_offset(sec_send) + bgq_weyl_section_offset(sec_recv);
	assert(offset_recv);
	return offset_recv;
}


static size_t bgq_offset_recv2send(size_t offset_recv) {
	bgq_weylfield_section sec_recv = bgq_sectionOfOffset(offset_recv);
	if ((sec_recv==sec_surface)||(sec_recv==sec_body))
		return offset_recv;
	//if (!COMM_T && ((sec_recv==sec_recv_tdown) || (sec_recv==sec_send_tup) || (sec_recv==sec_temp_tup)))
	//	return offset_recv;
	//if (!COMM_T && ((sec_recv==sec_recv_tup) || (sec_recv==sec_send_tdown)|| (sec_recv==sec_temp_tdown)))
	//	return offset_recv;
	bgq_weylfield_section sec_send = bgq_section_commbuftran(sec_recv, true);

	if (BGQ_UNVECTORIZE || !COMM_T) {
		if ((sec_recv==sec_recv_tdown) || (sec_recv==sec_send_tup) || (sec_recv==sec_temp_tup))
			sec_send = sec_temp_tup;
		if ((sec_recv==sec_recv_tup) || (sec_recv==sec_send_tdown)|| (sec_recv==sec_temp_tdown))
			return sec_temp_tdown;
	}

	size_t offset_send = offset_recv - bgq_weyl_section_offset(sec_recv) + bgq_weyl_section_offset(sec_send);
	assert(offset_send);
	return offset_send;
}


#if 0
// Find the offset in the weylfield where a specific value is expected for the HoppingMatrix kernel
// Used for 1st phase (distribute)
// i.e. body -> body,surface
//      surface -> body,surface or communication buffer to be sent to other node
static size_t bgq_weylfield_destoffsetForWeyl(bool isOdd, size_t ih_src, bgq_direction d_src) {
	bool isOdd_src = isOdd;
	bool isOdd_dst = !isOdd;
	size_t tv_src = bgq_halfvolume2tv(ih_src);
	size_t t1_src = bgq_halfvolume2t1(isOdd_src, ih_src);
	size_t t2_src = bgq_halfvolume2t2(isOdd_src, ih_src);
	size_t x_src = bgq_halfvolume2x(ih_src);
	size_t y_src = bgq_halfvolume2y(ih_src);
	size_t z_src = bgq_halfvolume2z(ih_src);
	//bool isSurface = bgq_physical2isSurface(isOdd_src, tv_src, x_src, y_src, z_src);
	//bgq_dimension dim = bgq_direction2dimension(d_src);

	size_t tv_dst = tv_src;
	size_t t1_dst = t1_src;
	size_t t2_dst = t2_src;
	size_t x_dst = x_src;
	size_t y_dst = y_src;
	size_t z_dst = z_src;
	bgq_direction d_dst = bgq_direction_revert(d_src);

	// Find the coordinate that consumes this value for its calculation
	// Change the coordinate to go into that direction
	switch (d_src) {
	case TDOWN:
		t1_dst = (t1_src + LOCAL_LT - 1) % LOCAL_LT;
		t2_dst = (t2_src + LOCAL_LT - 1) % LOCAL_LT;
		tv_dst = !bgq_physical2eo(isOdd_src, tv_src, x_src, y_src, z_src) ? (tv_src + PHYSICAL_LTV - 1) % PHYSICAL_LTV : tv_src;
		break;
	case TUP:
		t1_dst = (t1_src + 1) % LOCAL_LT;
		t2_dst = (t2_src + 1) % LOCAL_LT;
		tv_dst = !bgq_physical2eo(isOdd_src, tv_src, x_src, y_src, z_src) ? tv_src : (tv_src + 1) % PHYSICAL_LTV;
		break;
	case XDOWN:
		x_dst = (x_src + LOCAL_LX - 1) % LOCAL_LX;
		break;
	case XUP:
		x_dst = (x_src + 1) % LOCAL_LX;
		break;
	case YDOWN:
		y_dst = (y_src + LOCAL_LY - 1) % LOCAL_LY;
		break;
	case YUP:
		y_dst = (y_src + 1) % LOCAL_LY;
		break;
	case ZDOWN:
		z_dst = (z_src + LOCAL_LZ - 1) % LOCAL_LZ;
		break;
	case ZUP:
		z_dst = (z_src + 1) % LOCAL_LZ;
		break;
	}
	assert(bgq_local2tv(t1_dst,x_dst,y_dst,z_dst) == tv_dst);
	assert(bgq_local2tv(t2_dst,x_dst,y_dst,z_dst) == tv_dst);
	bgq_weylfield_section sec_read = bgq_HoppingMatrix_init_source_sectionof_physical(isOdd_dst, tv_dst, x_dst, y_dst, z_dst, d_dst); // dst is going to read the weyl
	assert(!bgq_section2isSend(sec_read));
	// Not that the section might logically belong to a different MPI node
	bgq_weylfield_section sec_write = bgq_section_commbuftran(sec_read, true); // src finds the correct location in preparation for dst
	assert(!bgq_section2isRecv(sec_write));
	// Find the offset at which the src node (can be this very same node) expects the datum
	size_t ih_dst = bgq_physical2halfvolume(tv_dst, x_dst, y_dst, z_dst);
	size_t offset_read = bgq_decode_offset(g_bgq_ih_dst2offset[isOdd][ih_dst].d[d_dst]);
	//assert(bgq_sectionOfOffset(offset_read) == sec_read);

	size_t offset_write;
	if ((sec_read == sec_body) || (sec_read == sec_surface)) {
		// No communication, do everything locally
		assert((sec_write==sec_body) || (sec_write==sec_surface));
		offset_write = offset_read;
	} else {
		// Communication with MPI buffers, translate recv- to send-buffer
		// i.e. dst expects the datum in the recv buffer, so we have to write it into the send buffer in order to be written there
		assert((sec_write!=sec_body) && (sec_write!=sec_surface));
		offset_write = bgq_offset_recv2send(offset_read);
	}
	//assert(offset_write);
	assert(0 <= offset_write && offset_write < bgq_weyl_section_offset(sec_end));
	assert(!bgq_section2isRecv(bgq_sectionOfOffset(offset_write)));
	return offset_write;
}
#endif

#if 0
size_t bgq_weyllayout_halfvolume2consecutiveoffset(bool isOdd_dst, size_t ih_dst, bgq_direction d_dst) {
	bool isSurface = bgq_halfvolume2isSurface(isOdd_dst, ih_dst);
	if (isSurface) {
		size_t is_dst = bgq_halfvolume2surface(isOdd_dst, ih_dst);
		return bgq_weyl_section_offset(sec_surface) + is_dst*sizeof(bgq_weylsite) + d_dst * sizeof(bgq_weyl_vec);
	} else {
		size_t ib_dst = bgq_halfvolume2body(isOdd_dst, ih_dst);
		return bgq_weyl_section_offset(sec_body) + ib_dst*sizeof(bgq_weylsite) + d_dst * sizeof(bgq_weyl_vec);
	}
}
#endif

bgq_direction bgq_offset2ddst(size_t offset) {
	bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
	switch (sec) {
	case sec_send_tup:
	case sec_recv_tdown:
	case sec_temp_tup:
		return TDOWN;
	case sec_send_tdown:
	case sec_recv_tup:
	case sec_temp_tdown:
		return TUP;
	case sec_send_xup:
	case sec_recv_xdown:
		return XDOWN;
	case sec_send_xdown:
	case sec_recv_xup:
		return XUP;
	case sec_send_yup:
	case sec_recv_ydown:
		return YDOWN;
	case sec_send_ydown:
	case sec_recv_yup:
		return YUP;
	case sec_send_zup:
	case sec_recv_zdown:
		return ZDOWN;
	case sec_send_zdown:
	case sec_recv_zup:
		return ZUP;
	case sec_surface:
	case sec_body: {
		size_t index = (offset - bgq_weyl_section_offset(sec_surface)) / sizeof(bgq_weyl_vec_double);
		return index % PHYSICAL_LD;
	}
	default:
		UNREACHABLE
	}
}


bgq_direction bgq_offset2dsrc(size_t offset) {
	bgq_direction d_dst = bgq_offset2ddst(offset);
	return bgq_direction_revert(d_dst);
}


size_t bgq_src2ih_dst(size_t t_src, size_t x_src, size_t y_src, size_t z_src, bgq_direction d_src) {
	size_t t_dst = t_src;
	size_t x_dst = x_src;
	size_t y_dst = y_src;
	size_t z_dst = z_src;

	switch (d_src) {
	case TDOWN:
		t_dst = (t_src + LOCAL_LT - 1) % LOCAL_LT;
		break;
	case TUP:
		t_dst = (t_src + 1) % LOCAL_LT;
		break;
	case XDOWN:
		x_dst = (x_src + LOCAL_LX - 1) % LOCAL_LX;
		break;
	case XUP:
		x_dst = (x_src + 1) % LOCAL_LX;
		break;
	case YDOWN:
		y_dst = (y_src + LOCAL_LY - 1) % LOCAL_LY;
		break;
	case YUP:
		y_dst = (y_src + 1) % LOCAL_LY;
		break;
	case ZDOWN:
		z_dst = (z_src + LOCAL_LZ - 1) % LOCAL_LZ;
		break;
	case ZUP:
		z_dst = (z_src + 1) % LOCAL_LZ;
		break;
	}

	size_t ih_dst = bgq_local2halfvolume(t_dst,x_dst,y_dst,z_dst);
	return ih_dst;
}


size_t bgq_src2k_dst(size_t t_src, size_t x_src, size_t y_src, size_t z_src, bgq_direction d_src) {
	size_t t_dst = t_src;
	size_t x_dst = x_src;
	size_t y_dst = y_src;
	size_t z_dst = z_src;

	switch (d_src) {
	case TDOWN:
		t_dst = (t_src + LOCAL_LT - 1) % LOCAL_LT;
		break;
	case TUP:
		t_dst = (t_src + 1) % LOCAL_LT;
		break;
	case XDOWN:
		x_dst = (x_src + LOCAL_LX - 1) % LOCAL_LX;
		break;
	case XUP:
		x_dst = (x_src + 1) % LOCAL_LX;
		break;
	case YDOWN:
		y_dst = (y_src + LOCAL_LY - 1) % LOCAL_LY;
		break;
	case YUP:
		y_dst = (y_src + 1) % LOCAL_LY;
		break;
	case ZDOWN:
		z_dst = (z_src + LOCAL_LZ - 1) % LOCAL_LZ;
		break;
	case ZUP:
		z_dst = (z_src + 1) % LOCAL_LZ;
		break;
	}

	size_t k_dst = bgq_local2k(t_dst, x_dst, y_dst, z_dst);
	return k_dst;
}


void bgq_indices_init() {
	if (g_bgq_indices_initialized)
		return;
	g_bgq_indices_initialized = true; // Take care for uses within this function itself

	g_comm_t = (g_nb_t_dn != g_proc_id);
	g_comm_x = (g_nb_x_dn != g_proc_id);
	g_comm_y = (g_nb_y_dn != g_proc_id);
	g_comm_z = (g_nb_z_dn != g_proc_id);
	assert(g_comm_t == (g_nb_t_up != g_proc_id));
	assert(g_comm_x == (g_nb_x_up != g_proc_id));
	assert(g_comm_y == (g_nb_y_up != g_proc_id));
	assert(g_comm_z == (g_nb_z_up != g_proc_id));
	//g_comm_t = true;
	//g_comm_x = true;
	//g_comm_y = true;
	//g_comm_z = true;
	g_bgq_dimension_isDistributed[DIM_T] = g_comm_t;
	g_bgq_dimension_isDistributed[DIM_X] = g_comm_x;
	g_bgq_dimension_isDistributed[DIM_Y] = g_comm_y;
	g_bgq_dimension_isDistributed[DIM_Z] = g_comm_z;
	g_bgq_dimension_hasHalo[DIM_T] = true;
	g_bgq_dimension_hasHalo[DIM_X] = g_comm_x;
	g_bgq_dimension_hasHalo[DIM_Y] = g_comm_y;
	g_bgq_dimension_hasHalo[DIM_Z] = g_comm_z;
	master_print("BGQ SPI/MPI communication enabled for: %s%s%s%s\n", COMM_T?"T,":"", COMM_X?"X,":"", COMM_Y?"Y,":"", COMM_Z?"Z":"");

	assert(PHYSICAL_LTV>=2);
	PHYSICAL_VOLUME = (PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY*PHYSICAL_LZ);
	if ((PHYSICAL_LTV <= 1) || (PHYSICAL_LX <= 2) || (PHYSICAL_LY <= 2) || (PHYSICAL_LZ <= 2)) {
		PHYSICAL_BODY = 0;
	} else {
		size_t body_tv = COMM_T ? (PHYSICAL_LTV - 1) : PHYSICAL_LTV;
		size_t body_x = COMM_X ? (PHYSICAL_LX - 2) : PHYSICAL_LX;
		size_t body_y = COMM_Y ? (PHYSICAL_LY - 2) : PHYSICAL_LY;
		size_t body_z = COMM_Z ? (PHYSICAL_LZ - 2) : PHYSICAL_LZ;
		PHYSICAL_BODY = body_tv * body_x * body_y * body_z;
	}
	PHYSICAL_SURFACE = PHYSICAL_VOLUME - PHYSICAL_BODY;

	// Setup array for index translation
	// It would be difficult to find an explicit expression for indices that are only on surface/body
	// Thus, we simply iterate over all locations and allocate locations in order
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		g_bgq_collapsed2halfvolume[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_collapsed2halfvolume[isOdd]));
		g_bgq_halfvolume2collapsed[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_halfvolume2collapsed[isOdd]));
#if 0
		g_bgq_index_surface2halfvolume[isOdd] = malloc(PHYSICAL_SURFACE * sizeof(*g_bgq_index_surface2halfvolume[isOdd]));
		g_bgq_index_body2halfvolume[isOdd] = malloc(PHYSICAL_BODY * sizeof(*g_bgq_index_body2halfvolume[isOdd]));
		g_bgq_index_halfvolume2surface[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_halfvolume2surface[isOdd]));
		g_bgq_index_halfvolume2body[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_halfvolume2body[isOdd]));
		g_bgq_index_halfvolume2surfacebody[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_halfvolume2surfacebody[isOdd]));
#endif
		size_t nextIndexSurface = 0;
		size_t nextIndexBody = 0;

		// Iteration order is effectively chosen here
		for (size_t ih = 0; ih < PHYSICAL_VOLUME; ih += 1) {
			size_t tv = bgq_halfvolume2tv(ih);
			size_t t1 = bgq_halfvolume2t1(isOdd, ih);
			size_t t2 = bgq_halfvolume2t2(isOdd, ih);
			size_t x = bgq_halfvolume2x(ih);
			size_t y = bgq_halfvolume2y(ih);
			size_t z = bgq_halfvolume2z(ih);
			assert(ih == bgq_local2halfvolume(t1,x,y,z));
			assert(ih == bgq_local2halfvolume(t2,x,y,z));
			assert(isOdd == bgq_local2isOdd(t1,x,y,z));
			assert(isOdd == bgq_local2isOdd(t2,x,y,z));

			bool isSurface = bgq_physical2isSurface(isOdd, tv, x, y, z);
			size_t ic;
			if (isSurface) {
				size_t is = nextIndexSurface;
				assert(0 <= is && is < PHYSICAL_SURFACE);
#if 0
				g_bgq_index_surface2halfvolume[isOdd][is] = ih;
				g_bgq_index_halfvolume2surfacebody[isOdd][ih] = is;
				g_bgq_index_halfvolume2surface[isOdd][ih] = is;
				g_bgq_index_halfvolume2body[isOdd][ih] = (size_t) -1;
#endif
				ic = bgq_surface2collapsed(is);
				nextIndexSurface += 1;
			} else {
				size_t ib = nextIndexBody;
				assert(0 <= ib && ib < PHYSICAL_BODY);
#if 0
				g_bgq_index_body2halfvolume[isOdd][ib] = ih;
				g_bgq_index_halfvolume2surfacebody[isOdd][ih] = ib;
				g_bgq_index_halfvolume2surface[isOdd][ih] = (size_t) -1;
				g_bgq_index_halfvolume2body[isOdd][ih] = ib;
#endif
				ic = bgq_body2collapsed(ib);
				nextIndexBody += 1;
			}
			g_bgq_collapsed2halfvolume[isOdd][ic] = ih;
			g_bgq_halfvolume2collapsed[isOdd][ih] = ic;
		}
	}

#if 0
#ifndef NDEBUG
	// Check consistency
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		for (size_t is2 = 0; is2 < PHYSICAL_SURFACE; is2 += 1) {
			size_t ih = bgq_surface2halfvolume(isOdd, is2);
			assert(bgq_halfvolume2isSurface(isOdd,ih));

			assert(is2 == bgq_halfvolume2surface(isOdd,ih));
			assert(ih == bgq_surface2halfvolume(isOdd,is2));
		}
	}
#endif
#endif

	// Setup mapping of weyl to some memory offset
	// (where to read a datum for hoppingmatrix)
	// Note: although we arrange data ordered as in ih_dst, the field contains data ih_src
	assert(bgq_weyl_section_offset(sec_end) % sizeof(bgq_weyl_vec_double) == 0);
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		size_t indices = bgq_weyl_section_offset(sec_end) / sizeof(bgq_weyl_vec_double);
		g_bgq_index2collapsed[isOdd] = malloc(indices * sizeof(*g_bgq_index2collapsed[isOdd]));
		//g_bgq_collapsed2indexrecv[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_collapsed2indexsend[isOdd]));
		g_bgq_collapsed2indexsend[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_collapsed2indexsend[isOdd]));

#if 0
		g_bgq_index2ih_dst[isOdd] = malloc(indices * sizeof(*g_bgq_index2ih_dst[isOdd]));
		g_bgq_index2d_dst[isOdd] = malloc(indices * sizeof(*g_bgq_index2d_dst[isOdd]));
		g_bgq_ih_dst2offset[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_ih_dst2offset[isOdd]));
		g_bgq_is_dst2offset[isOdd] = malloc(PHYSICAL_SURFACE * sizeof(*g_bgq_is_dst2offset[isOdd]));
		g_bgq_ib_dst2offset[isOdd] = malloc(PHYSICAL_BODY * sizeof(*g_bgq_ib_dst2offset[isOdd]));
#endif
	}
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		bool isOdd_src = !isOdd;
		bool isOdd_dst = isOdd;

		size_t nextoffset[sec_end];
		for (bgq_weylfield_section sec = 0; sec < sec_end; sec += 1) {
			nextoffset[sec] = bgq_weyl_section_offset(sec);
		}
		//if (!COMM_T) {
		//	// In this case, send==recv
		//	nextoffset[sec_recv_tup] = nextoffset[sec_send_tdown];
		//	nextoffset[sec_recv_tdown] = nextoffset[sec_send_tup];
		//}

		for (ucoord ih_dst = 0; ih_dst < PHYSICAL_VOLUME; ih_dst += 1) {
			ucoord ic_dst = bgq_halfvolume2collapsed(isOdd_dst, ih_dst);
			size_t tv_dst = bgq_halfvolume2tv(ih_dst);
			size_t x_dst = bgq_halfvolume2x(ih_dst);
			size_t y_dst = bgq_halfvolume2y(ih_dst);
			size_t z_dst = bgq_halfvolume2z(ih_dst);
			bool isSurface_dst = bgq_collapsed2isSurface(ic_dst);
			bgq_weylfield_section mainsec_dst = isSurface_dst ? sec_surface : sec_body;

			for (size_t d_dst = 0; d_dst < PHYSICAL_LD; d_dst += 1) {
				bgq_direction d_src = bgq_direction_revert(d_dst);
				bgq_dimension dim = bgq_direction2dimension(d_src);
				ucoord ic_src = bgq_collapsed_dst2src(isOdd_dst, ic_dst, d_dst);
				ucoord ih_src = bgq_collapsed2halfvolume(isOdd_src, ic_src);
				size_t tv_src = bgq_halfvolume2tv(ih_src);
				size_t x_src = bgq_halfvolume2x(ih_src);
				size_t y_src = bgq_halfvolume2y(ih_src);
				size_t z_src = bgq_halfvolume2z(ih_src);
				bool isSurface_src = bgq_collapsed2isSurface(ic_src);

				//bgq_weylfield_section mainsec_src = isSurface_src ? sec_surface : sec_body;
				bgq_weylfield_section sec_write = bgq_HoppingMatrix_init_source_sectionof_physical(isOdd_src, tv_src, x_src, y_src, z_src, d_src);

				// Reserve some offset in surface/body to ensure consecutive layout
				size_t offset_main = nextoffset[mainsec_dst];
				nextoffset[mainsec_dst] += sizeof(bgq_weyl_vec_double);
				assert(bgq_collapsed2consecutiveoffset(ic_dst, d_dst) == offset_main);

				ucoord index_main = bgq_offset2index(offset_main);
				if (ic_dst == 60) {
					int a = 0;
				}
				g_bgq_index2collapsed[isOdd_dst][index_main] = ic_dst;
				assert(bgq_sectionOfOffset(offset_main) == mainsec_dst);

				if (sec_write != mainsec_dst) {
					// If in one of the comm send buffers, also reserve some space there
					assert((sec_write!=sec_surface) && (sec_write!=sec_body));

					size_t offset_write = nextoffset[sec_write];
					nextoffset[sec_write] += sizeof(bgq_weyl_vec_double);

					ucoord index_write = bgq_offset2index(offset_write);
					g_bgq_collapsed2indexsend[isOdd_dst][ic_src/*!!!*/].d[d_dst] = index_write;
					if (ic_dst == 60) {
						int a = 0;
					}
					g_bgq_index2collapsed[isOdd_dst][index_write] = ic_dst;

#if 0
					// Sender part (this is on the inter-node sending the data to this node)
					size_t offset_sender = bgq_offset_recv2send(thisOffset);
					ucoord index_send = bgq_offset2index(offset_sender);
					//assert(bgq_sectionOfOffset(offset_sender) == bgq_section_commbuftran(sec,true));
					g_bgq_index2collapsed[isOdd_dst][index_send] = ic_dst;
					g_bgq_collapsed2indexsend[isOdd_dst][ic_src].d[d_src] = index_send;
					//g_bgq_index2ih_dst[isOdd][offset_sender / sizeof(bgq_weyl_vec)] = ih_dst;
					//g_bgq_index2d_dst[isOdd][offset_sender / sizeof(bgq_weyl_vec)] = d_dst;
#endif
				} else {
					g_bgq_collapsed2indexsend[isOdd_dst][ic_src/*!!!*/].d[d_dst] = index_main;
					//g_bgq_collapsed2indexrecv[isOdd_dst][ic_dst].d[d_dst] = index_main;
				}

#if 0
				//assert(bgq_sectionOfOffset(thisOffset) == sec);
				g_bgq_ih_dst2offset[isOdd][ih_dst].d[d_dst] = bgq_encode_offset(thisOffset);
				if (isSurface) {
					size_t is_dst = bgq_halfvolume2surface(isOdd_dst, ih_dst);
					g_bgq_is_dst2offset[isOdd][is_dst].d[d_dst] = bgq_encode_offset(thisOffset);
				} else {
					size_t ib_dst = bgq_halfvolume2body(isOdd_dst, ih_dst);
					g_bgq_ib_dst2offset[isOdd][ib_dst].d[d_dst] = bgq_encode_offset(thisOffset);
				}
#endif
			}
		}
	}

#if 0
#ifndef NDEBUG
	// Check consistency
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		bool isOdd_dst = !isOdd;
		bool isOdd_src = isOdd;

		for (size_t is_dst = 0; is_dst < PHYSICAL_SURFACE; is_dst += 1) {
			for (size_t d_dst = 0; d_dst < PHYSICAL_LD; d_dst += 1) {
				size_t offset_is = g_bgq_is_dst2offset[isOdd][is_dst].d[d_dst];
				//assert(offset_is);
				// For valgrind to check initialization

				size_t ih_dst = bgq_surface2halfvolume(isOdd_dst, is_dst);
				size_t offset_ih = g_bgq_ih_dst2offset[isOdd][ih_dst].d[d_dst];
				//assert(offset_ih);

				assert(offset_is == offset_ih);
			}
		}
	}
#endif
#endif

#if 0
	// Determine dest locations (where to write a weyl before communication)
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		g_bgq_ihsrc2offsetwrite[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_ihsrc2offsetwrite[isOdd]));
		g_bgq_issrc2offsetwrite[isOdd] = malloc(PHYSICAL_SURFACE * sizeof(*g_bgq_issrc2offsetwrite[isOdd]));
		g_bgq_ibsrc2offsetwrite[isOdd] = malloc(PHYSICAL_BODY * sizeof(*g_bgq_ibsrc2offsetwrite[isOdd]));

		size_t indices = bgq_weyl_section_offset(sec_end) / sizeof(bgq_weyl_vec);
		g_bgq_index2ih_src[isOdd] = malloc(indices * sizeof(*g_bgq_index2ih_src[isOdd]));
		g_bgq_index2d_src[isOdd] = malloc(indices * sizeof(*g_bgq_index2d_src[isOdd]));
	}
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		bool isOdd_src = isOdd;
		bool isOdd_dst = !isOdd;
		for (size_t ih_src = 0; ih_src < PHYSICAL_VOLUME; ih_src += 1) {
			for (size_t d_src = 0; d_src < PHYSICAL_LD; d_src += 1) {
				size_t offset_send = bgq_weylfield_destoffsetForWeyl(isOdd, ih_src, d_src);
				size_t index_send = offset_send / sizeof(bgq_weyl_vec);

				assert(!bgq_section2isRecv(bgq_sectionOfOffset(offset_send)));

				size_t encodedOffset = bgq_encode_offset(offset_send);
				g_bgq_ihsrc2offsetwrite[isOdd][ih_src].d[d_src] = encodedOffset;
				bool isSurface = bgq_halfvolume2isSurface(isOdd_src, ih_src);
				if (isSurface) {
					size_t is_src = bgq_halfvolume2surface(isOdd_src, ih_src);
					g_bgq_issrc2offsetwrite[isOdd][is_src].d[d_src] = encodedOffset;
				} else {
					size_t ib_src = bgq_halfvolume2body(isOdd, ih_src);
					g_bgq_ibsrc2offsetwrite[isOdd][ib_src].d[d_src] = encodedOffset;
				}

				// Reverse lookup
				// 3 different offset calculations: send,recv,consecutive; may overlap some regions
				size_t index_write = offset_send / sizeof(bgq_weyl_vec);
				g_bgq_index2ih_src[isOdd][index_write] = ih_src;
				g_bgq_index2d_src[isOdd][index_write] = d_src;
				assert(d_src == bgq_offset2dsrc(offset_send));

				size_t offset_receiver = bgq_offset_send2recv(offset_send);
				assert(!bgq_section2isRecv(offset_receiver));
				assert(bgq_offset_recv2send(offset_receiver) == offset_send);
				size_t index_receiver = offset_receiver / sizeof(bgq_weyl_vec);
				g_bgq_index2ih_src[isOdd][index_receiver] = ih_src;
				g_bgq_index2d_src[isOdd][index_receiver] = d_src;
				assert(d_src == bgq_offset2dsrc(offset_receiver));

				size_t tv_src = bgq_halfvolume2tv(ih_src);
				size_t t1_src = bgq_halfvolume2t1(isOdd_dst, ih_src);
				size_t t2_src = bgq_halfvolume2t2(isOdd_dst, ih_src);
				size_t x_src = bgq_halfvolume2x(ih_src);
				size_t y_src = bgq_halfvolume2y(ih_src);
				size_t z_src = bgq_halfvolume2z(ih_src);
				bgq_direction d_dst = bgq_direction_revert(d_src);
				size_t ih_dst = bgq_src2ih_dst(t1_src, x_src, y_src, z_src, d_src);
				assert(ih_dst == bgq_src2ih_dst(t2_src, x_src, y_src, z_src, d_src));
				size_t offset_consecutive = bgq_weyllayout_halfvolume2consecutiveoffset(isOdd_dst,ih_dst,d_dst);
				size_t index_consecutive = offset_consecutive / sizeof(bgq_weyl_vec);
				g_bgq_index2ih_src[isOdd][index_consecutive] = ih_src;
				g_bgq_index2d_src[isOdd][index_consecutive] = d_src;
			}
		}
	}

#ifndef NDEBUG
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		for (size_t offset = bgq_weyl_section_offset(0); offset < bgq_weyl_section_offset(sec_end); offset += sizeof(bgq_weyl_vec)) {
			size_t index = bgq_offset2index(offset);
			bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
			if (bgq_section2isRecv(sec)) {
				size_t offsend = bgq_offset_recv2send(offset);
			}


			size_t ih_dst = g_bgq_index2ih_dst[isOdd][index];
			assert(0 <= ih_dst && ih_dst < PHYSICAL_VOLUME);
			bgq_direction d_dst = g_bgq_index2d_dst[isOdd][index];
			assert(d_dst == bgq_offset2ddst(offset));

			size_t ih_src = g_bgq_index2ih_src[isOdd][index];
			assert(0 <= ih_src && ih_src < PHYSICAL_VOLUME);
			bgq_direction d_src = g_bgq_index2d_src[isOdd][index];
			assert(d_src == bgq_offset2dsrc(offset));

			assert(d_dst == bgq_direction_revert(d_src));
			//assert(bgq_src2ih_dst());
		}
	}
#endif
#endif
}



