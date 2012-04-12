#pragma once

//nomenclature: 
//-glb is relative to the global grid
//-loc to the local one
int glb_size[4],glb_vol,glb_spat_vol,glb_volh;
int loc_size[4],loc_vol,loc_spat_vol,loc_volh;
int *dir_of_bord;
//-lx is lexicografic
coords *glb_coord_of_loclx;
coords *loc_coord_of_loclx;
int *glblx_of_loclx;
int *glblx_of_bordlx;
int *loclx_of_bordlx;
int *glblx_of_edgelx;
int nissa_lx_geom_inited;
//-eo is even-odd
int *loclx_parity;
int *loceo_of_loclx;
int *loclx_of_loceo[2];
coords *loceo_neighup[2];
coords *loceo_neighdw[2];
int nissa_eo_geom_inited;

//neighbours of local volume + borders
coords *loclx_neighdw,*loclx_neighup;
coords *loclx_neigh[2];

//stopping criterions for multimass inverter
#define sc_standard 0
#define sc_unilevel 1
#define sc_weighted_norm2 2
#define sc_weighted_norm_inf 3

#define numb_known_stopping_criterion 4
char list_known_stopping_criterion[4][1024]={"standard","unilevel","weighted_norm2","weighted_norm_inf"};

//error handler
MPI_Errhandler mpi_error_handler;

//basic mpi types
MPI_Datatype MPI_FLOAT_128;
MPI_Datatype MPI_SU3;
MPI_Datatype MPI_QUAD_SU3;
MPI_Datatype MPI_COLOR;
MPI_Datatype MPI_SPINCOLOR;
MPI_Datatype MPI_SPINCOLOR_128;
MPI_Datatype MPI_REDSPINCOLOR;

//float 128 summ
MPI_Op MPI_FLOAT_128_SUM;

//timings
double tot_nissa_time=0;
double tot_nissa_comm_time=0;

//verbosity
const double nissa_default_verbosity=2;
double nissa_verbosity=2;

//size of the border along the 4 dir,types for sending
int start_lx_bord_send_up[4],start_lx_bord_rece_up[4];
int start_lx_bord_send_dw[4],start_lx_bord_rece_dw[4];
int bord_dir_vol[4],bord_offset[4];
int bord_vol,bord_volh;
MPI_Datatype MPI_LX_SU3_BORDS_SEND[4],MPI_LX_SU3_BORDS_RECE[4];
MPI_Datatype MPI_LX_QUAD_SU3_BORDS_SEND[4],MPI_LX_QUAD_SU3_BORDS_RECE[4];
MPI_Datatype MPI_LX_SPINCOLOR_BORDS_SEND[4],MPI_LX_SPINCOLOR_BORDS_RECE[4];
MPI_Datatype MPI_LX_SPINCOLOR_128_BORDS_SEND[4],MPI_LX_SPINCOLOR_128_BORDS_RECE[4];
int start_eo_bord_send_up[4],start_eo_bord_rece_up[4];//,start_ev_bord_send_up,start_od_bord_send_up;
int start_eo_bord_send_dw[4],start_eo_bord_rece_dw[4];//,start_ev_bord_send_dw,start_od_bord_send_dw;
MPI_Datatype MPI_EO_QUAD_SU3_BORDS_SEND_TXY[4],MPI_EO_QUAD_SU3_BORDS_RECE[4];
MPI_Datatype MPI_EV_QUAD_SU3_BORDS_SEND_Z[2],MPI_OD_QUAD_SU3_BORDS_SEND_Z[2];
MPI_Datatype MPI_EO_COLOR_BORDS_SEND_TXY[4],MPI_EO_COLOR_BORDS_RECE[4];
MPI_Datatype MPI_EV_COLOR_BORDS_SEND_Z[2],MPI_OD_COLOR_BORDS_SEND_Z[2];
MPI_Datatype MPI_EO_SPINCOLOR_BORDS_SEND_TXY[4],MPI_EO_SPINCOLOR_BORDS_RECE[4];
MPI_Datatype MPI_EV_SPINCOLOR_BORDS_SEND_Z[2],MPI_OD_SPINCOLOR_BORDS_SEND_Z[2];

int bord_offset_eo[2][8]; //eo, 8 dirs

//size of the edges along the 6 directions
int edge_dir_vol[6],edge_offset[6];
int edge_vol,edge_volh;
int edge_numb[4][4]={{-1,0,1,2},{0,-1,3,4},{1,3,-1,5},{2,4,5,-1}};
MPI_Datatype MPI_LX_SU3_EDGES_SEND[6],MPI_LX_SU3_EDGES_RECE[6];
MPI_Datatype MPI_LX_QUAD_SU3_EDGES_SEND[6],MPI_LX_QUAD_SU3_EDGES_RECE[6];
MPI_Datatype MPI_EO_QUAD_SU3_EDGES_SEND[96],MPI_EO_QUAD_SU3_EDGES_RECE[6];

int rank,rank_tot,cart_rank;
int nparal_dir;
coords paral_dir;
coords nrank_dir;
coords rank_coord;
coords rank_neigh[2],rank_neighdw,rank_neighup;
coords plan_rank,line_rank;

int little_endian;

#define nreals_per_color 6
#define nreals_per_spincolor 24
#define nreals_per_quad_su3 72

MPI_Comm cart_comm;
MPI_Comm plan_comm[4];
MPI_Comm line_comm[4];

//vectors
int nissa_max_required_memory;
int nissa_required_memory;
void *main_nissa_arr;
nissa_vect main_nissa_vect;
nissa_vect *last_nissa_vect;

//random generator stuff
rnd_gen glb_rnd_gen;
int nissa_glb_rnd_gen_inited;
rnd_gen *loc_rnd_gen;
int nissa_loc_rnd_gen_inited;

as2t smunu_entr[4];   //these are the sigma matrices entries
int smunu_pos[4][6];  //and positions

//the real amd imaginary unit
complex ONE;
complex I;

//The base of the 16 gamma matrixes and the two rotators
dirac_matr base_gamma[16];
dirac_matr Pplus,Pminus;
const char gtag[16][3]={"S0","V1","V2","V3","V0","P5","A1","A2","A3","A0","T1","T2","T3","B1","B2","B3"};
