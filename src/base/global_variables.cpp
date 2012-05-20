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
EXTERN int *dir_of_bord;
//-lx is lexicografic
EXTERN coords *glb_coord_of_loclx;
EXTERN coords *loc_coord_of_loclx;
EXTERN int *glblx_of_loclx;
EXTERN int *glblx_of_bordlx;
EXTERN int *loclx_of_bordlx;
EXTERN int *glblx_of_edgelx;
EXTERN int nissa_lx_geom_inited;
//-eo is even-odd
EXTERN int *loclx_parity;
EXTERN int *loceo_of_loclx;
EXTERN int *loclx_of_loceo[2];
EXTERN coords *loceo_neighup[2];
EXTERN coords *loceo_neighdw[2];
EXTERN int nissa_eo_geom_inited;

//neighbours of local volume + borders
EXTERN coords *loclx_neighdw,*loclx_neighup;
EXTERN coords *loclx_neigh[2];

//basic mpi types
EXTERN MPI_Datatype MPI_FLOAT_128;
EXTERN MPI_Datatype MPI_SU3;
EXTERN MPI_Datatype MPI_QUAD_SU3;
EXTERN MPI_Datatype MPI_COLOR;
EXTERN MPI_Datatype MPI_SPIN;
EXTERN MPI_Datatype MPI_SPINCOLOR;
EXTERN MPI_Datatype MPI_SPINCOLOR_128;
EXTERN MPI_Datatype MPI_REDSPINCOLOR;

//float 128 summ
EXTERN MPI_Op MPI_FLOAT_128_SUM;

//timings
EXTERN double tot_nissa_time;
EXTERN double tot_nissa_comm_time;

//verbosity
EXTERN int verb_call;
EXTERN int nissa_verbosity;

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

EXTERN int bord_offset_eo[2][8]; //eo, 8 dirs

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
EXTERN int rank,rank_tot,cart_rank;
EXTERN int nparal_dir;
EXTERN coords paral_dir;
EXTERN coords nrank_dir;
EXTERN coords rank_coord;
EXTERN coords rank_neigh[2],rank_neighdw,rank_neighup;
EXTERN coords plan_rank,line_rank,line_coord_rank;
EXTERN int nissa_grid_inited;

//endianess
EXTERN int little_endian;

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

//random generator stuff
EXTERN rnd_gen glb_rnd_gen;
EXTERN int nissa_glb_rnd_gen_inited;
EXTERN rnd_gen *loc_rnd_gen;
EXTERN int nissa_loc_rnd_gen_inited;

EXTERN as2t smunu_entr[4];   //these are the sigma matrices entries
EXTERN int smunu_pos[4][6];  //and positions

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
