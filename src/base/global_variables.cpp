#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../new_types/new_types_definitions.h"

#ifdef ONLY_INSTANTIATION
 #define EXTERN extern
#else
 #define EXTERN
#endif

//nomenclature: 
//-glb is relative to the global grid
//-loc to the local one
EXTERN int glb_size[4],glb_vol,glb_spat_vol,glb_volh;
EXTERN int loc_size[4],loc_vol,loc_spat_vol,loc_volh;
EXTERN int bulk_vol,non_bw_surf_vol,non_fw_surf_vol;
EXTERN int surf_vol,bw_surf_vol,fw_surf_vol;
EXTERN double glb_vol2,loc_vol2;
//-lx is lexicografic
EXTERN coords *glb_coord_of_loclx;
EXTERN coords *loc_coord_of_loclx;
EXTERN int *glblx_of_loclx;
EXTERN int *glblx_of_bordlx;
EXTERN int *loclx_of_bordlx;
EXTERN int *surflx_of_bordlx;
EXTERN int *glblx_of_edgelx;
EXTERN int *loclx_of_bulklx;
EXTERN int *loclx_of_surflx;
EXTERN int *loclx_of_non_bw_surflx;
EXTERN int *loclx_of_non_fw_surflx;
EXTERN int *loclx_of_bw_surflx;
EXTERN int *loclx_of_fw_surflx;
EXTERN int nissa_lx_geom_inited;
//-eo is even-odd
EXTERN int *loclx_parity;
EXTERN int *loceo_of_loclx;
EXTERN int *loclx_of_loceo[2];
EXTERN int *surfeo_of_bordeo[2];
EXTERN coords *loceo_neighup[2];
EXTERN coords *loceo_neighdw[2];
EXTERN int nissa_eo_geom_inited;
EXTERN int nissa_use_eo_geom;

//neighbours of local volume + borders
EXTERN coords *loclx_neighdw,*loclx_neighup;
EXTERN coords *loclx_neigh[2];

//basic mpi types
EXTERN MPI_Datatype MPI_FLOAT_128;
EXTERN MPI_Datatype MPI_SU3;
EXTERN MPI_Datatype MPI_QUAD_SU3;
EXTERN MPI_Datatype MPI_COLOR;
EXTERN MPI_Datatype MPI_SPIN;
EXTERN MPI_Datatype MPI_SPINSPIN;
EXTERN MPI_Datatype MPI_SPINCOLOR;
EXTERN MPI_Datatype MPI_SPINCOLOR_128;
EXTERN MPI_Datatype MPI_REDSPINCOLOR;

//float 128 summ
EXTERN MPI_Op MPI_FLOAT_128_SUM;

//timings
EXTERN int ncgm_inv,ncg_inv;
EXTERN double cgm_inv_over_time,cg_inv_over_time;
EXTERN double tot_nissa_time;
EXTERN double tot_nissa_comm_time;

//verbosity
EXTERN int verb_call;
EXTERN int nissa_verbosity;
EXTERN int nissa_warn_if_not_disallocated;
EXTERN int nissa_warn_if_not_communicated;

//128 bit precision
EXTERN int nissa_use_128_bit_precision;

//size of the border along the 4 dir,types for sending
EXTERN int start_lx_bord_send_up[4],start_lx_bord_rece_up[4];
EXTERN int start_lx_bord_send_dw[4],start_lx_bord_rece_dw[4];
EXTERN int bord_dir_vol[4],bord_offset[4];
EXTERN int bord_vol,bord_volh;
EXTERN MPI_Datatype MPI_LX_SU3_BORDS_SEND[4],MPI_LX_SU3_BORDS_RECE[4];
EXTERN MPI_Datatype MPI_LX_QUAD_SU3_BORDS_SEND[4],MPI_LX_QUAD_SU3_BORDS_RECE[4];
EXTERN MPI_Datatype MPI_LX_SPIN_BORDS_SEND[4],MPI_LX_SPIN_BORDS_RECE[4];
EXTERN MPI_Datatype MPI_LX_COLOR_BORDS_SEND[4],MPI_LX_COLOR_BORDS_RECE[4];
EXTERN MPI_Datatype MPI_LX_SPINSPIN_BORDS_SEND[4],MPI_LX_SPINSPIN_BORDS_RECE[4];
EXTERN MPI_Datatype MPI_LX_SPINCOLOR_BORDS_SEND[4],MPI_LX_SPINCOLOR_BORDS_RECE[4];
EXTERN MPI_Datatype MPI_LX_SPINCOLOR_128_BORDS_SEND[4],MPI_LX_SPINCOLOR_128_BORDS_RECE[4];
EXTERN int start_eo_bord_send_up[4],start_eo_bord_rece_up[4];
EXTERN int start_eo_bord_send_dw[4],start_eo_bord_rece_dw[4];
EXTERN MPI_Datatype MPI_EO_QUAD_SU3_BORDS_SEND_TXY[4],MPI_EO_QUAD_SU3_BORDS_RECE[4];
EXTERN MPI_Datatype MPI_EV_QUAD_SU3_BORDS_SEND_Z[2],MPI_OD_QUAD_SU3_BORDS_SEND_Z[2];
EXTERN MPI_Datatype MPI_EO_COLOR_BORDS_SEND_TXY[4],MPI_EO_COLOR_BORDS_RECE[4];
EXTERN MPI_Datatype MPI_EV_COLOR_BORDS_SEND_Z[2],MPI_OD_COLOR_BORDS_SEND_Z[2];
EXTERN MPI_Datatype MPI_EO_SPIN_BORDS_SEND_TXY[4],MPI_EO_SPIN_BORDS_RECE[4];
EXTERN MPI_Datatype MPI_EV_SPIN_BORDS_SEND_Z[2],MPI_OD_SPIN_BORDS_SEND_Z[2];
EXTERN MPI_Datatype MPI_EO_SPINCOLOR_BORDS_SEND_TXY[4],MPI_EO_SPINCOLOR_BORDS_RECE[4];
EXTERN MPI_Datatype MPI_EV_SPINCOLOR_BORDS_SEND_Z[2],MPI_OD_SPINCOLOR_BORDS_SEND_Z[2];
EXTERN MPI_Datatype MPI_EO_SPINCOLOR_128_BORDS_SEND_TXY[4],MPI_EO_SPINCOLOR_128_BORDS_RECE[4];
EXTERN MPI_Datatype MPI_EV_SPINCOLOR_128_BORDS_SEND_Z[2],MPI_OD_SPINCOLOR_128_BORDS_SEND_Z[2];

EXTERN int bord_offset_eo[2][8]; //eo, 8 dirs

//use async communication
EXTERN int nissa_use_async_communications;

//size of the edges along the 6 directions
EXTERN int edge_dir_vol[6],edge_offset[6];
EXTERN int edge_vol,edge_volh;
EXTERN int edge_numb[4][4]
#ifndef ONLY_INSTANTIATION
={{-1,0,1,2},{0,-1,3,4},{1,3,-1,5},{2,4,5,-1}}
#endif
;
EXTERN MPI_Datatype MPI_LX_SU3_EDGES_SEND[6],MPI_LX_SU3_EDGES_RECE[6];
EXTERN MPI_Datatype MPI_LX_QUAD_SU3_EDGES_SEND[6],MPI_LX_QUAD_SU3_EDGES_RECE[6];
EXTERN MPI_Datatype MPI_EO_QUAD_SU3_EDGES_SEND[96],MPI_EO_QUAD_SU3_EDGES_RECE[6];

//ranks
EXTERN int rank,nissa_nranks,cart_rank;
EXTERN int nparal_dir;
EXTERN coords paral_dir;
EXTERN coords nrank_dir;
EXTERN coords rank_coord;
EXTERN coords rank_neigh[2],rank_neighdw,rank_neighup;
EXTERN coords plan_rank,line_rank,line_coord_rank;
EXTERN int nissa_grid_inited;

//thread
#ifdef THREAD_DEBUG
EXTERN int glb_barr_line;
EXTERN char *glb_barr_file;
#endif
#ifndef ONLY_INSTANTIATION
int thread_pool_locked=true,nthreads=1;
#else
EXTERN int thread_pool_locked,nthreads;
#endif
EXTERN double *glb_double_reduction_buf;
EXTERN float_128 *glb_float_128_reduction_buf;

EXTERN void(*threaded_function_ptr)();

//endianess
EXTERN int little_endian;

//global input file handle
EXTERN FILE *input_global;

//volume, plan and line communicator
EXTERN MPI_Comm cart_comm;
EXTERN MPI_Comm plan_comm[4];
EXTERN MPI_Comm line_comm[4];

//vectors
EXTERN int nissa_max_required_memory;
EXTERN int nissa_required_memory;
EXTERN void *main_nissa_arr;
EXTERN nissa_vect main_nissa_vect;
EXTERN nissa_vect *last_nissa_vect;
EXTERN void *return_nissa_malloc_ptr;

//random generator stuff
EXTERN rnd_gen glb_rnd_gen;
EXTERN int nissa_glb_rnd_gen_inited;
EXTERN rnd_gen *loc_rnd_gen;
EXTERN int nissa_loc_rnd_gen_inited;
EXTERN enum rnd_t nissa_rnd_type_map[5]
#ifndef ONLY_INSTANTIATION
={RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_Z2,RND_Z2,RND_Z4}
#endif
;
EXTERN as2t smunu_entr[4];   //these are the sigma matrices entries
EXTERN int smunu_pos[4][6];  //and positions

//perpendicular dir
#ifdef ONLY_INSTANTIATION
EXTERN int perp_dir[4][3],perp2_dir[4][3][2],perp3_dir[4][3][2];
#else
int perp_dir[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
int perp2_dir[4][3][2]={{{2,3},{1,3},{1,2}},{{2,3},{0,3},{2,3}},{{1,3},{0,3},{0,1}},{{1,2},{0,2},{0,1}}};
int perp3_dir[4][3][2]={{{3,2},{3,1},{2,1}},{{3,2},{3,0},{3,2}},{{3,1},{3,0},{1,0}},{{2,1},{2,0},{1,0}}};
#endif

//the real amd imaginary unit
EXTERN complex ONE,I;

//The base of the 16 gamma matrixes and the two rotators
EXTERN dirac_matr base_gamma[16];
EXTERN dirac_matr Pplus,Pminus;
EXTERN char gtag[16][3]
#ifndef ONLY_INSTANTIATION
={"S0","V1","V2","V3","V0","P5","A1","A2","A3","A0","T1","T2","T3","B1","B2","B3"}
#endif
;
EXTERN int nissa_map_mu[4]
#ifndef ONLY_INSTANTIATION
={4,1,2,3};
#endif
;
EXTERN spinspin nissa_opg[4],nissa_omg[4];

EXTERN int su3_sub_gr_indices[3][2]
#ifndef ONLY_INSTANTIATION
={{0,1},{1,2},{0,2}};
#endif
;

/////////////////////////////////////////// buffered comm ///////////////////////////////////

EXTERN int nbuffered_comm_allocated;
EXTERN int buffered_comm_in_prog;

//buffers
EXTERN size_t nissa_buff_size;
EXTERN char *nissa_recv_buf,*nissa_send_buf;

//buffered communicators
EXTERN buffered_comm_t buffered_lx_spincolor_comm;
EXTERN buffered_comm_t buffered_lx_halfspincolor_comm;
EXTERN buffered_comm_t buffered_lx_colorspinspin_comm;
EXTERN buffered_comm_t buffered_lx_su3spinspin_comm;
EXTERN buffered_comm_t buffered_lx_quad_su3_comm;
EXTERN buffered_comm_t buffered_eo_color_comm;
EXTERN buffered_comm_t buffered_eo_quad_su3_comm;

/////////////////////////////////////////// BGQ specifics ///////////////////////////////////

#ifdef BGQ

//indices of remapping and output hopping matrix
EXTERN int *bgqlx_of_loclx,*loclx_of_bgqlx;
EXTERN bi_halfspincolor *bgq_hopping_matrix_output_data;
EXTERN bi_halfspincolor **bgq_hopping_matrix_output_pointer;
EXTERN bi_halfspincolor *bgq_hopping_matrix_output_T_buffer;
EXTERN int bgqlx_t_vbord_vol,bgq_vsurf_vol;

#endif

/////////////////////////////////////////// SPI specifics ///////////////////////////////////

#ifdef SPI

#include <spi/include/kernel/MU.h>
#include <spi/include/mu/InjFifo.h>
#include <spi/include/mu/GIBarrier.h>

//flag to remember if spi has been initialized
EXTERN int nissa_spi_inited;

//spi rank coordinates
EXTERN coords_5D spi_rank_coord;

//neighbours in the 4 dirs
EXTERN MUHWI_Destination_t spi_neigh[2][4];

//spi fifo and counters for bytes
EXTERN uint64_t *spi_fifo[8],spi_desc_count[8];
EXTERN MUSPI_InjFifoSubGroup_t spi_fifo_sg_ptr;

//spi barrier
EXTERN MUSPI_GIBarrier_t spi_barrier;

//bats
EXTERN MUSPI_BaseAddressTableSubGroup_t spi_bat_gr;
EXTERN uint32_t spi_bat_id[2];

//counter
EXTERN volatile uint64_t spi_recv_counter;

//physical address
EXTERN uint64_t spi_send_buf_phys_addr;

#endif

#undef EXTERN
