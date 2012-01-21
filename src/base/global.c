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

//stopping criterions for multimass inverter
#define sc_standard 0
#define sc_unilevel 1
#define sc_weighted_norm2 2
#define sc_weighted_norm_inf 3

#define numb_known_stopping_criterion 4
char list_known_stopping_criterion[4][1024]={"standard","unilevel","weighted_norm2","weighted_norm_inf"};

//basic mpi types
MPI_Datatype MPI_SU3;
MPI_Datatype MPI_QUAD_SU3;
MPI_Datatype MPI_COLOR;
MPI_Datatype MPI_SPINCOLOR;
MPI_Datatype MPI_REDSPINCOLOR;

//size of the border along the 4 dir,types for sending
int start_lx_bord_send_up[4],start_lx_bord_rece_up[4];
int start_lx_bord_send_dw[4],start_lx_bord_rece_dw[4];
int bord_dir_vol[4],bord_offset[4];
int loc_bord,loc_bordh;
MPI_Datatype MPI_LX_SU3_BORD_SEND[4],MPI_LX_SU3_BORD_RECE[4];
MPI_Datatype MPI_LX_GAUGE_BORD_SEND[4],MPI_LX_GAUGE_BORD_RECE[4];
MPI_Datatype MPI_LX_SPINCOLOR_BORD_SEND[4],MPI_LX_SPINCOLOR_BORD_RECE[4];
int start_eo_bord_send_up[4],start_eo_bord_rece_up[4];
int start_eo_bord_send_dw[4],start_eo_bord_rece_dw[4];
MPI_Datatype MPI_EO_GAUGE_BORD_SEND[4],MPI_EO_GAUGE_BORD_RECE[4];
MPI_Datatype MPI_EO_COLOR_BORD_SEND[4],MPI_EO_COLOR_BORD_RECE[4];
//MPI_Datatype MPI_LXREDSPINCOLOR_BORD[4];

int bord_offset_eo[2][8]; //eo, 8 dirs

//size of the edges along the 6 directions
int edge_dir_vol[6],edge_offset[6];
int loc_edge,loc_edgeh;
int edge_numb[4][4]={{-1,0,1,2},{0,-1,3,4},{1,3,-1,5},{2,4,5,-1}};
MPI_Datatype MPI_LX_GAUGE_EDGE_SEND[6],MPI_LX_GAUGE_EDGE_RECE[6];
MPI_Datatype MPI_EO_GAUGE_EDGE_SEND[6],MPI_EO_GAUGE_EDGE_RECE[6];

int paral_dir[4],nparal_dir;
int nrank_dir[4];
int rank_coord[4];
int rank_neighdw[4],rank_neighup[4];
int rank,rank_tot,cart_rank;
int plan_rank[4],line_rank[4];

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
