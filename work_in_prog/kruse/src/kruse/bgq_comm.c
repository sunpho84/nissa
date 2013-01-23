/*
 * bgq_comm.c
 *
 *  Created on: Oct 25, 2012
 *      Author: meinersbur
 */

#define BGQ_COMM_C_
#include "bgq_comm.h"

#include "bgq_field.h"
#include "bgq_dispatch.h"
#include "bgq_workers.h"

#ifdef SPI
#include "../DirectPut.h"
#include <upci/upc_atomic.h>
#endif
#include "../global.h"

#include <mpi.h>




static MPI_Request g_bgq_request_recv_double[PHYSICAL_LD];
static MPI_Request g_bgq_request_send_double[PHYSICAL_LD];

static MPI_Request g_bgq_request_recv_float[PHYSICAL_LD];
static MPI_Request g_bgq_request_send_float[PHYSICAL_LD];



static bgq_dimension bgq_commdim2dimension(ucoord commdim) {
	assert(0 <= commdim && commdim < PHYSICAL_LD);
	if (COMM_T) {
		if (commdim == 0)
			return DIM_T;
		commdim -= 1;
	}
	if (COMM_X) {
		if (commdim == 0)
			return DIM_X;
		commdim -= 1;
	}
	if (COMM_Y) {
		if (commdim == 0)
			return DIM_Y;
		commdim -= 1;
	}
	assert(COMM_Z);
	assert(commdim==0);
	return DIM_Z;
}


static ucoord bgq_direction2commdir(bgq_direction d) {
	ucoord result = 0;
	if (COMM_T) {
		if (d==TUP)
			return result;
		result+=1;
		if (d==TDOWN)
			return result;
		result+=1;
	}
	if (COMM_X) {
		if(d==XUP)
			return result;
		result+=1;
		if (d==XDOWN)
			return result;
		result+=1;
	}
	if (COMM_Y) {
		if(d==YUP)
			return result;
		result+=1;
		if (d==YDOWN)
			return result;
		result+=1;
	}
	if (COMM_Z) {
		if(d==ZUP)
			return result;
		result+=1;
		if (d==ZDOWN)
			return result;
		result+=1;
	}
	UNREACHABLE
	return -1;
}


static ucoord bgq_commdir2direction(ucoord commdir) {
	if (COMM_T) {
		if (commdir==0)
			return TUP;
		commdir-=1;
		if (commdir==0)
			return TDOWN;
		commdir-=1;
	}
	if (COMM_X) {
		if (commdir==0)
			return XUP;
		commdir-=1;
		if (commdir==0)
			return XDOWN;
		commdir-=1;
	}
	if (COMM_Y) {
		if (commdir==0)
			return YUP;
		commdir-=1;
		if (commdir==0)
			return YDOWN;
		commdir-=1;
	}
	if (COMM_Z) {
		if (commdir==0)
			return ZUP;
		commdir-=1;
		if (commdir==0)
			return ZDOWN;
		commdir-=1;
	}
	UNREACHABLE
	return -1;
}


static int bgq_direction2rank(bgq_direction d) {
	switch (d) {
	case TUP:
		return g_nb_t_up;
	case TDOWN:
		return g_nb_t_dn;
	case XUP:
		return g_nb_x_up;
	case XDOWN:
		return g_nb_x_dn;
	case YUP:
		return g_nb_y_up;
	case YDOWN:
		return g_nb_y_dn;
	case ZUP:
		return g_nb_z_up;
	case ZDOWN:
		return g_nb_z_dn;
	default:
		UNREACHABLE
	}
	return -1;
}


static MPI_Request g_bgq_request_recv_double[PHYSICAL_LD];
static MPI_Request g_bgq_request_send_double[PHYSICAL_LD];

static MPI_Request g_bgq_request_recv_float[PHYSICAL_LD];
static MPI_Request g_bgq_request_send_float[PHYSICAL_LD];


#ifdef SPI
static unsigned mySpirank;
static inline unsigned bgq_abcde2spirank(Personality_t *pers, uint8_t a, uint8_t b,uint8_t c,uint8_t d,uint8_t e) {
	assert(pers);
	torus_t tdims = {
		pers->Network_Config.Anodes,
		pers->Network_Config.Bnodes,
		pers->Network_Config.Cnodes,
		pers->Network_Config.Dnodes,
		pers->Network_Config.Enodes
	};
	torus_t dims = {
		pers->Network_Config.Anodes,
		pers->Network_Config.Bnodes,
		pers->Network_Config.Cnodes,
		pers->Network_Config.Dnodes,
		pers->Network_Config.Enodes
	};

	unsigned numNodes = tdims.a * tdims.b * tdims.c * tdims.d * tdims.e;
	unsigned result = ((((a)*dims.b + b)*dims.c + c)*dims.d + d)*dims.e + e;
	assert(result < numNodes);
	return result;
}


static void setup_destinations(Personality_t *pers) {
	torus_t tcoords = {
		pers->Network_Config.Acoord,
		pers->Network_Config.Bcoord,
		pers->Network_Config.Ccoord,
		pers->Network_Config.Dcoord,
		pers->Network_Config.Ecoord
	};

	torus_t tdims = {
		pers->Network_Config.Anodes,
		pers->Network_Config.Bnodes,
		pers->Network_Config.Cnodes,
		pers->Network_Config.Dnodes,
		pers->Network_Config.Enodes
	};

	//numNodes = tdims.a * tdims.b * tdims.c * tdims.d * tdims.e;
	mySpirank = bgq_abcde2spirank(pers, tcoords.a, tcoords.b, tcoords.c, tcoords.e, tcoords.d);
	//myRank = g_proc_id;
	//myRank = mySpirank;

	if (Kernel_ProcessCount() > 1) {
		printf("ERROR: this test only works with ppn==1\n");
		abort();
	}

	for (ucoord cd = 0; cd < COMMDIR_COUNT; cd += 1) {
		bgq_direction d_src = bgq_commdir2direction(cd);
		bgq_direction d_dst = bgq_direction_revert(d_src);

		MPI_Status status;
		int mpirank_to = bgq_direction2rank(d_src);
		int mpirank_from = bgq_direction2rank(d_dst);
		torus_t nb = {-1,-1,-1,-1,-1};
		MPI_CHECK(MPI_Sendrecv(&tcoords, sizeof(tcoords), MPI_BYTE, mpirank_to, 0, &nb, sizeof(nb), MPI_BYTE, mpirank_from, 0, g_cart_grid, &status));
		assert(get_MPI_count(&status) == sizeof(tcoords));

		ucoord nbrank = bgq_abcde2spirank(pers, nb.a, nb.b, nb.c, nb.d, nb.e);
		nb2dest[cd].tcoord = nb; // Just for reference
		MUSPI_SetUpDestination(&nb2dest[cd].dest, nb.a, nb.b, nb.c, nb.d, nb.e);
		nb2dest[cd].hintsABCD = 0;
		nb2dest[cd].hintsE = 0;
		//printf("node %d: %d(%d,%d,%d,%d,%d)-%d->%d(%d,%d,%d,%d,%d)\n", g_proc_id, mySpirank, tcoords.a, tcoords.b, tcoords.c, tcoords.e, tcoords.d, cd, nbrank, nb.a, nb.b, nb.c, nb.d, nb.e);

		//nb2dest[COMMDIR_COUNT + cd] = nb2dest[cd];
	}
}
#endif



static void bgq_comm_test(bool nospi, bool sloppy) {
	// test communication

	for (ucoord cd = 0; cd < COMMDIR_COUNT; cd+=1) {
		bgq_direction d = bgq_commdir2direction(cd);
		bgq_weylfield_section sec = bgq_direction2section(d, true);
		ucoord index_end = bgq_section_size(sec) / sizeof(bgq_weyl_vec_double);
		for (size_t i = 0; i < index_end; i += 1) {
			for (ucoord v = 0; v < 2; v += 1){
				for (ucoord c = 0; c < 3; c += 1){
					for (ucoord k = 0; k < 2; k += 1){
						if (sloppy)
							g_bgq_sec_send_float[d][i].s[v][c][k] = g_cart_id + g_proc_id * _Complex_I;
						else
							g_bgq_sec_send_double[d][i].s[v][c][k] = g_cart_id + g_proc_id * _Complex_I;

						g_bgq_sec_recv_double[d][i].s[v][c][k] = -1;
						g_bgq_sec_recv_float[d][i].s[v][c][k] = -1;
					}
				}
			}
		}

	}



	bgq_comm_recv(nospi, sloppy, NULL);
	bgq_comm_send();
	bgq_comm_wait();

	for (ucoord cd = 0; cd < COMMDIR_COUNT; cd+=1) {
		bgq_direction d = bgq_commdir2direction(cd);
		assert(bgq_direction_isDistributed(d));

		bgq_weylfield_section sec = bgq_direction2section(d, false);
		ucoord index_end = bgq_section_size(sec) / sizeof(bgq_weyl_vec_double);
		int rank_neighbor = bgq_direction2rank(d);
		for (size_t i = 0; i < index_end; i += 1) {
			for (ucoord v = 0; v < 2; v += 1) {
				for (ucoord c = 0; c < 3; c += 1) {
					for (ucoord k = 0; k < 2; k += 1) {
						complexdouble val = sloppy ?  g_bgq_sec_recv_float[d][i].s[v][c][k] : g_bgq_sec_recv_double[d][i].s[v][c][k];
						if (creal(val) != rank_neighbor) {
							printf("Node %d: Exchange doesn't work for dir %d: %d != %f (proc=%f) at point %zu\n", g_proc_id, d, rank_neighbor, creal(val), cimag(val), i);
							abort();
						}
					}
				}
			}
		}
	}

	master_print("%s Communication sloppy=%d tested successfully\n", nospi ? "MPI" : "MU SPI", sloppy);
}


static bool g_bgq_comm_common_initialized = false;
static void bgq_comm_common_init(void) {
	if (g_bgq_comm_common_initialized)
		return;
	g_bgq_comm_common_initialized = true;
	bgq_indices_init();

	if (PHYSICAL_SURFACE > 0) {
		size_t commbufsize = bgq_weyl_section_offset(sec_comm_end) - bgq_weyl_section_offset(sec_comm);
		uint8_t *buf = (uint8_t*)malloc_aligned(commbufsize /*+ commbufsize    some memory is wasted here, but no additional work to be done for alignment*/, BGQ_ALIGNMENT_L2);
		uint8_t *bufend = buf + commbufsize;
		g_bgq_sec_comm = buf;
		g_bgq_sec_comm_float = buf + commbufsize;
		for (bgq_direction d = 0; d < PHYSICAL_LD; d += 1) {
			bgq_dimension dim = bgq_direction2dimension(d);

			bgq_weylfield_section sec_send = bgq_direction2section(d, true);
			bgq_weylfield_section sec_recv = bgq_direction2section(d, false);

			g_bgq_sec_send_double[d] = (bgq_weyl_vec_double*)(buf + bgq_weyl_section_offset(sec_send) - bgq_weyl_section_offset(sec_comm));
			assert((uint8_t*)g_bgq_sec_send_double[d] <= bufend);
			g_bgq_sec_send_float[d] = (bgq_weyl_vec_float*)((uint8_t*)g_bgq_sec_send_double[d] + 0);

			g_bgq_sec_recv_double[d] = (bgq_weyl_vec_double*)(buf + bgq_weyl_section_offset(sec_recv) - bgq_weyl_section_offset(sec_comm));
			assert((uint8_t*)g_bgq_sec_recv_double[d] <= bufend);
			g_bgq_sec_recv_float[d] = (bgq_weyl_vec_float*)((uint8_t*)g_bgq_sec_recv_double[d] + 0);

			assert((uintptr_t)g_bgq_sec_send_double[d] % BGQ_ALIGNMENT_L2 == 0);
			assert((uintptr_t)g_bgq_sec_recv_double[d] % BGQ_ALIGNMENT_L2 == 0);
		}
	}

	if (BGQ_UNVECTORIZE || !COMM_T) {
		g_bgq_sec_temp_tup_double = malloc_aligned(LOCAL_HALO_T/PHYSICAL_LP *sizeof(bgq_weyl_vec_double), BGQ_ALIGNMENT_L2);
		g_bgq_sec_temp_tup_float = (bgq_weyl_vec_float*)g_bgq_sec_temp_tup_double;
		g_bgq_sec_temp_tdown_double = malloc_aligned(LOCAL_HALO_T/PHYSICAL_LP *sizeof(bgq_weyl_vec_double), BGQ_ALIGNMENT_L2);
		g_bgq_sec_temp_tdown_float = (bgq_weyl_vec_float*)g_bgq_sec_temp_tdown_double;
		//g_bgq_sec_vrecv_tup = malloc_aligned(PHYSICAL_HALO_T,BGQ_ALIGNMENT_L2);
		//g_bgq_sec_vrecv_tdown = malloc_aligned(PHYSICAL_HALO_T,BGQ_ALIGNMENT_L2);
	}
}


static bool g_bgq_comm_mpi_initialized = false;
void bgq_comm_mpi_init(void) {
	if (g_bgq_comm_mpi_initialized)
		return;
	g_bgq_comm_mpi_initialized = true;
	bgq_comm_common_init();
	if (PHYSICAL_SURFACE == 0)
		return;

	for (ucoord cd = 0; cd < COMMDIR_COUNT; cd+=1) {
		bgq_direction d_src = bgq_commdir2direction(cd);
		assert(bgq_direction_isDistributed(d_src));

		//master_print("d_src=%d\n",d_src);
		bgq_direction d_dst = bgq_direction_revert(d_src);
		//master_print("d_dst=%d\n",d_dst);
		bgq_dimension dim = bgq_direction2dimension(d_src);
		//master_print("dim=%d\n",dim);

		ucoord commdir_src = cd;
		//master_print("commdir_src=%d\n",commdir_src);
		bgq_weylfield_section sec_recv = bgq_direction2section(d_src, false);
		size_t secsize = bgq_weyl_section_offset(sec_recv+1) - bgq_weyl_section_offset(sec_recv);
		assert(secsize > 0);
		MPI_CHECK(MPI_Recv_init(g_bgq_sec_recv_double[d_src], secsize / sizeof(complex_double), MPI_DOUBLE_COMPLEX, bgq_direction2rank(d_src), d_dst, g_cart_grid, &g_bgq_request_recv_double[cd]));
		MPI_CHECK(MPI_Recv_init(g_bgq_sec_recv_float[d_src], secsize / sizeof(complex_double), MPI_COMPLEX, bgq_direction2rank(d_src), d_dst, g_cart_grid, &g_bgq_request_recv_float[cd]));
		//master_print("MPI_CHECK(MPI_Recv_init(%zu, %zu, %zu, %zu, %zu, %zu, %zu))\n", (size_t)(g_bgq_sec_recv[d_src]), secsize / sizeof(double), (size_t)(MPI_DOUBLE), (size_t)bgq_direction2rank(d_src), (size_t)d_dst, (size_t)g_cart_grid, (size_t)(&g_bgq_request_recv[commdir_src]));

		ucoord commdir_dst = bgq_direction2commdir(d_dst);
		bgq_weylfield_section sec_send = bgq_direction2section(d_dst, true);
		assert(secsize == bgq_weyl_section_offset(sec_send+1) - bgq_weyl_section_offset(sec_send));
		MPI_CHECK(MPI_Send_init(g_bgq_sec_send_double[d_dst], secsize / sizeof(complex_double), MPI_DOUBLE_COMPLEX, bgq_direction2rank(d_dst), d_dst, g_cart_grid, &g_bgq_request_send_double[cd]));
		MPI_CHECK(MPI_Send_init(g_bgq_sec_send_float[d_dst], secsize / sizeof(complex_double), MPI_COMPLEX, bgq_direction2rank(d_dst), d_dst, g_cart_grid, &g_bgq_request_send_float[cd]));
		//master_print("MPI_CHECK(MPI_Send_init(%zu, %zu, %zu, %zu, %zu, %zu, %zu))\n", g_bgq_sec_send[d_dst], secsize / sizeof(double), MPI_DOUBLE, bgq_direction2rank(d_dst), d_dst, g_cart_grid, &g_bgq_request_send[commdir_dst]);
	}

	bgq_comm_test(true, false);
	bgq_comm_test(true, true);
}


static uint64_t totalMessageSize_float;

static bool g_bgq_comm_spi_initialized = false;
void bgq_comm_spi_init(void) {
#ifdef SPI
	if (g_bgq_comm_spi_initialized)
		return;
	g_bgq_comm_spi_initialized = true;
	bgq_comm_common_init();
	if (PHYSICAL_SURFACE == 0)
		return;

	size_t messageSizes[2*COMMDIR_COUNT];
	size_t roffsets[2*COMMDIR_COUNT];
	size_t soffsets[2*COMMDIR_COUNT];

	size_t messageSizes_float[COMMDIR_COUNT];
	//size_t roffsets_flaot[COMMDIR_COUNT];
	//size_t soffsets_float[COMMDIR_COUNT];

	// here comes the SPI initialization
	int spi_num_dirs = COMMDIR_COUNT;
	totalMessageSize = 0;
	totalMessageSize_float = 0;
	size_t bufsize_double = bgq_weyl_section_offset(sec_send_end) - bgq_weyl_section_offset(sec_send_begin);
	for (ucoord cd = 0; cd < COMMDIR_COUNT; cd+=1) {
		bgq_direction d_src = bgq_commdir2direction(cd);
		bgq_direction d_dst = bgq_direction_revert(d_src);
		bgq_dimension dim = bgq_direction2dimension(d_src);
		assert(bgq_direction_isDistributed(d_src));

		ucoord commdir_src = bgq_direction2commdir(d_src);
		bgq_weylfield_section sec_recv = bgq_direction2section(d_src, false);
		size_t secsize = bgq_section_size(sec_recv);
		ucoord commdir = commdir_src;

		ucoord commdir_dst = bgq_direction2commdir(d_dst);
		bgq_weylfield_section sec_send = bgq_direction2section(d_dst, true);
		assert(secsize == bgq_section_size(sec_send));


		messageSizes[commdir] = secsize;
		soffsets[commdir] = bgq_weyl_section_offset(sec_send) - bgq_weyl_section_offset(sec_send_begin);
		assert(soffsets[commdir] % BGQ_ALIGNMENT_L2 == 0);
		assert(soffsets[commdir] + messageSizes[commdir] <= bgq_weyl_section_offset(sec_send_end) - bgq_weyl_section_offset(sec_send_begin));
		roffsets[cd] = bgq_weyl_section_offset(sec_recv) - bgq_weyl_section_offset(sec_recv_begin);
		assert(roffsets[cd] % BGQ_ALIGNMENT_L2 == 0);
		assert((roffsets[cd] + secsize) <= (bgq_weyl_section_offset(sec_recv_end) - bgq_weyl_section_offset(sec_recv_begin)));
		totalMessageSize += secsize;


		size_t secsize_float = secsize/2;
		messageSizes_float[commdir] = secsize_float;
		messageSizes[COMMDIR_COUNT+commdir] = secsize_float;
		//soffsets_float[commdir] = bufsize_double + soffsets[commdir];
		soffsets[COMMDIR_COUNT+commdir] = bufsize_double + soffsets[commdir];
		//roffsets_commdir[cd] = bufsize_double + roffsets[cd];
		roffsets[COMMDIR_COUNT+cd] = bufsize_double + roffsets[cd];
		totalMessageSize_float += secsize_float;

		//master_print("SPI %llu: d=%llu msize=%zu soffset=%zu d_dst=%llu roffset=%zu\n", cd, d_src, messageSizes[commdir], soffsets[commdir], d_dst, roffsets[cd]);
	}
	assert(totalMessageSize == bgq_weyl_section_offset(sec_recv_end) - bgq_weyl_section_offset(sec_recv_begin));
	assert(totalMessageSize_float == totalMessageSize/2);






	do_dynamic = 0; // Use static routing (since only neighbor-to-neighbor communication)


	//do_dynamic = 0; // Use static routing (since only neighbor-to-neighbor communication)
	//char *argv[] = {"exec"};
	//main2(1, argv);
	//master_print("Main2 done\n");

	Personality_t pers;
	int rc = 0;
	// get the CNK personality
	Kernel_GetPersonality(&pers, sizeof(pers));
	setup_destinations(&pers);

	// adjust the SPI pointers to the send and receive buffers
	SPIrecvBuffers = (char*)g_bgq_sec_recv_double[0];
	assert((uintptr_t)SPIrecvBuffers % BGQ_ALIGNMENT_L2 == 0);
	SPIsendBuffers = (char*)g_bgq_sec_send_double[0];
	assert((uintptr_t)SPIsendBuffers % BGQ_ALIGNMENT_L2 == 0);

	// Setup the FIFO handles
	rc = msg_InjFifoInit(&injFifoHandle,
			 0,                      /* startingSubgroupId */
			 0,                      /* startingFifoId     */
			 COMMDIR_COUNT,          /* numFifos   */
			 INJ_MEMORY_FIFO_SIZE+1, /* fifoSize */
			 NULL                    /* Use default attributes */
			 );
	if(rc != 0) {
		fprintf(stderr, "msg_InjFifoInit failed with rc=%d\n", rc);
		exit(1);
	}

	// Set up base address table for reception counter and buffer
	uint64_t bufferSize = bgq_weyl_section_offset(sec_recv_end) - bgq_weyl_section_offset(sec_recv_begin);
	assert(bufferSize == totalMessageSize);
	setup_mregions_bats_counters(bufferSize);

	// Create descriptors
	// Injection Direct Put Descriptor, one for each neighbor
	SPIDescriptors = (MUHWI_Descriptor_t*)(((uint64_t)SPIDescriptorsMemory+64)&~(64-1));
	create_descriptors(SPIDescriptors, messageSizes, soffsets, roffsets, COMMDIR_COUNT);

	SPIDescriptors32 = (MUHWI_Descriptor_t*)(((uint64_t)SPIDescriptorsMemory32+64)&~(64-1));
	create_descriptors(SPIDescriptors32, messageSizes_float, soffsets, roffsets, COMMDIR_COUNT);

    // Initialize the barrier, resetting the hardware.
    rc = MUSPI_GIBarrierInit(&GIBarrier, 0 /*comm world class route */);
    if (rc) {
    	printf("MUSPI_GIBarrierInit returned rc=%d\n",rc);
    	abort();
    }

#if 0
    for (size_t cd = 0; cd < COMMDIR_COUNT; cd+=1) {
    	bgq_direction d = bgq_commdir2direction(cd);
    	bgq_weylfield_section sec_send = bgq_direction2section(d,true);
    	size_t secsize = bgq_section_size(sec_send);

    	SPIDescriptors[cd].Message_Length = messageSizes[cd];
    	SPIDescriptors[cd].Pa_Payload     = sendBufPAddr + soffsets[cd];
    	MUSPI_SetRecPutOffset(&SPIDescriptors[cd], roffsets[cd]);
    }
#endif


	bgq_comm_test(false, false);
	bgq_comm_test(false, true);
#endif
}


static bool g_bgq_comm_mpiActive = false;
static bool g_bgq_comm_mpiSloppy;
static bgq_weylfield_controlblock *g_bgq_comm_mpiTarget;
#ifdef SPI
static bool g_bgq_comm_spiActive = false;
static bool g_bgq_comm_spiSloppy;
static bgq_weylfield_controlblock *g_bgq_comm_spiTarget;
#endif

//TODO: inline?
void bgq_comm_recv(bool nospi, bool sloppy, bgq_weylfield_controlblock *targetfield) {
	if (PHYSICAL_SURFACE == 0)
		return;
	assert(omp_get_thread_num()==0);
	//master_print("Comm Receiving... nospi=%d sloppy=%d\n", nospi, sloppy);
	bgq_comm_wait(); // Ensure that previous communication has finished

#ifdef SPI
	if (!nospi) {
	    // reset the recv counter
	    recvCounter = sloppy ? totalMessageSize_float : totalMessageSize;
	    g_bgq_comm_spiActive = true;
	    g_bgq_comm_spiSloppy = sloppy;
	    g_bgq_comm_spiTarget = targetfield;
		return;
	}
#endif

	if (sloppy)
		MPI_CHECK(MPI_Startall(COMMDIR_COUNT, g_bgq_request_recv_float));
	else
		MPI_CHECK(MPI_Startall(COMMDIR_COUNT, g_bgq_request_recv_double));
	g_bgq_comm_mpiActive = true;
    g_bgq_comm_mpiSloppy = sloppy;
    g_bgq_comm_mpiTarget = targetfield;
}


void bgq_comm_send(void) {
	if (PHYSICAL_SURFACE == 0)
		return;
	assert(omp_get_thread_num()==0);
	//master_print("Comm Sending... nospi=%d sloppy=%d\n", nospi, sloppy);

#ifdef SPI
	assert(g_bgq_comm_spiActive || g_bgq_comm_mpiActive);
	if (g_bgq_comm_spiActive) {
		// make sure everybody has reset recvCounter
		global_barrier();

		for (size_t cd = 0; cd < COMMDIR_COUNT; cd += 1) {
			descCount[cd] = msg_InjFifoInject(injFifoHandle, cd, g_bgq_comm_spiSloppy ? &SPIDescriptors32[cd] : &SPIDescriptors[cd]);
			if (descCount[cd] == -1) {
				printf("msg_InjFifoInject failed, most likely because there is no room in the fifo\n");
				abort();
			}
		}
	}
#endif

	if (g_bgq_comm_mpiActive) {
		if (g_bgq_comm_mpiSloppy)
			MPI_CHECK(MPI_Startall(COMMDIR_COUNT, g_bgq_request_send_float));
		else
			MPI_CHECK(MPI_Startall(COMMDIR_COUNT, g_bgq_request_send_double));
	}
}


static void bgq_comm_datamove_field(bgq_weylfield_controlblock *targetfield) {
	assert(PHYSICAL_SURFACE > 0);
	if (!targetfield)
		return; // bgq_comm_test

	// 5. Move received to correct location, so the communication buffers can be reused
	if (targetfield->pendingDatamove) {
		bgq_master_sync();
		static bgq_work_datamove work;
		work.spinorfield = targetfield;
		work.opts = targetfield->hmflags;
		assert(targetfield->has_weyllayout_float + targetfield->has_weyllayout_double == 1);
		if (targetfield->has_weyllayout_float)
			bgq_master_call(&bgq_HoppingMatrix_worker_datamove_float, &work);
		if (targetfield->has_weyllayout_double)
			bgq_master_call(&bgq_HoppingMatrix_worker_datamove_double, &work);
		targetfield->pendingDatamove = false;
	}
}


void bgq_comm_wait(void) {
	if (PHYSICAL_SURFACE == 0)
		return;
	assert(omp_get_thread_num()==0);
	//master_print("Comm Waiting... nospi=%d sloppy=%d\n", nospi, sloppy);

#ifdef SPI
	if (g_bgq_comm_spiActive) {
		uint64_t startTime = 0;
		uint64_t expectedBytes = g_bgq_comm_spiSloppy ? totalMessageSize_float : totalMessageSize;

		// Wait for all data is received
		 //printf("node %d: %llu bytes to be received\n", g_proc_id, expectedBytes);
		 while(recvCounter > 0) {
			 // Check range of pending bytes to receive
			 assert(recvCounter <= expectedBytes);

#if 0
			 if (GetTimeBase() - startTime >= 1600) {
				 printf("node %d: %llu bytes left\n", g_proc_id, recvCounter);
				 startTime = GetTimeBase();
			 }
#endif
		 }
		 //printf("node %d: All data received\n", g_proc_id);

		 // Wait for all data sent
		while (true) {
			size_t sendDone = 0;
			for (size_t j = 0; j < COMMDIR_COUNT; j += 1) {
				sendDone += msg_InjFifoCheckCompletion(injFifoHandle, j, descCount[j]);
			}
			if (sendDone == COMMDIR_COUNT)
				break;
		}

		_bgq_msync();  // Ensure data is available to all cores.
		bgq_comm_datamove_field(g_bgq_comm_spiTarget);
		g_bgq_comm_spiTarget = NULL;
		g_bgq_comm_spiActive = false;
	}
#endif

	if (g_bgq_comm_mpiActive) {
		MPI_Status recv_status[COMMDIR_COUNT];
		if (g_bgq_comm_mpiSloppy)
			MPI_CHECK(MPI_Waitall(COMMDIR_COUNT, g_bgq_request_recv_float, recv_status));
		else
			MPI_CHECK(MPI_Waitall(COMMDIR_COUNT, g_bgq_request_recv_double, recv_status));

		MPI_Status send_status[COMMDIR_COUNT];
		if (g_bgq_comm_mpiSloppy)
			MPI_CHECK(MPI_Waitall(COMMDIR_COUNT, g_bgq_request_send_float, send_status));
		else
			MPI_CHECK(MPI_Waitall(COMMDIR_COUNT, g_bgq_request_send_double, send_status));

#ifndef NDEBUG
		for (ucoord commdir = 0; commdir < COMMDIR_COUNT; commdir += 1) {
			bgq_direction d = bgq_commdir2direction(commdir);
			bgq_weylfield_section sec = bgq_direction2section(d, false);
			size_t size = bgq_weyl_section_offset(sec+1) - bgq_weyl_section_offset(sec);
			if (g_bgq_comm_mpiSloppy)
				size /= 2;
			assert(get_MPI_count(&recv_status[commdir]) == size);
		}
#endif

		bgq_comm_datamove_field(g_bgq_comm_mpiTarget);
		g_bgq_comm_mpiTarget = NULL;
		g_bgq_comm_mpiActive = false;
	}
}


