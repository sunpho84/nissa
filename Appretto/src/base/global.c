#pragma once

#ifdef APPRETTO
#define EXTERN 
#else
#define EXTERN extern
#endif

//nomenclature: 
//-glb is relative to the global grid
//-loc to the local one
EXTERN int glb_size[4],glb_vol;
EXTERN int loc_size[4],loc_vol;
EXTERN int *dir_of_bord;
//-lx is lexicografic
EXTERN int **glb_coord_of_loclx;
EXTERN int **loc_coord_of_loclx;
EXTERN int *glblx_of_loclx;
EXTERN int *glblx_of_bordlx;
EXTERN int *loclx_of_bordlx;
EXTERN int *glblx_of_edgelx;
//-eo is even-odd
EXTERN int *loclx_parity;
EXTERN int *loceo_of_loclx;
EXTERN int *loclx_of_loce;
EXTERN int *loclx_of_loco;
EXTERN int **loce_neighup;
EXTERN int **loce_neighdw;
EXTERN int **loco_neighup;
EXTERN int **loco_neighdw;
EXTERN int loc_volr;
EXTERN int appretto_eo_geom_init;

//neighbours of local volume + borders
EXTERN int **loclx_neighdw,**loclx_neighup;

//stopping criterions for multimass inverter
#define sc_standard 0
#define sc_unilevel 1
#define sc_weighted_norm2 2
#define sc_weighted_norm_inf 3

#define numb_known_stopping_criterion 4
#ifndef APPRETTO
char list_known_stopping_criterion[4][1024]={"standard","unilevel","weighted_norm2","weighted_norm_inf"};
#endif

//basic mpi types
EXTERN MPI_Datatype MPI_SU3;
EXTERN MPI_Datatype MPI_QUAD_SU3;
EXTERN MPI_Datatype MPI_SPINCOLOR;
EXTERN MPI_Datatype MPI_REDSPINCOLOR;

//size of the border along the 4 dir,types for sending
EXTERN int start_lx_bord_send_up[4],start_lx_bord_rece_up[4];
EXTERN int start_lx_bord_send_dw[4],start_lx_bord_rece_dw[4];
EXTERN int bord_dir_vol[4],bord_offset[4];
EXTERN int loc_bord;
EXTERN MPI_Datatype MPI_SU3_BORD_SEND[4],MPI_SU3_BORD_RECE[4];
EXTERN MPI_Datatype MPI_GAUGE_BORD_SEND[4],MPI_GAUGE_BORD_RECE[4];
EXTERN MPI_Datatype MPI_LXSPINCOLOR_BORD_SEND[4],MPI_LXSPINCOLOR_BORD_RECE[4];
EXTERN MPI_Datatype MPI_LXREDSPINCOLOR_BORD[4];

EXTERN int bord_offset_eo[2][8]; //eo, 8 dirs

//size of the edges along the 6 directions
EXTERN int edge_dir_vol[6],edge_offset[6];
EXTERN int loc_edge;
#ifdef APPRETTO
int edge_numb[4][4]={{-1,0,1,2},{0,-1,3,4},{1,3,-1,5},{2,4,5,-1}};
#endif
EXTERN MPI_Datatype MPI_GAUGE_EDGE_SEND[6],MPI_GAUGE_EDGE_RECE[6];

EXTERN int paral_dir[4],nparal_dir;
EXTERN int nproc_dir[4];
EXTERN int proc_coord[4];
EXTERN int rank_neighdw[4],rank_neighup[4];
EXTERN int rank,rank_tot,cart_rank;

EXTERN int big_endian;

#define nreals_per_spincolor 24
#define nreals_per_quad_su3 72

EXTERN MPI_Comm cart_comm;
EXTERN MPI_Comm plan_comm;
EXTERN MPI_Comm time_comm;

//vectors
EXTERN int appretto_max_required_memory;
EXTERN int appretto_required_memory;
EXTERN void *main_appretto_arr;
EXTERN appretto_vect main_appretto_vect;
EXTERN appretto_vect *last_appretto_vect;

//random generator stuff
#define ran2_ntab 32
EXTERN int *ran2_idum,*ran2_idum2,**ran2_iv,*ran2_iy;
EXTERN int random_initialized;

EXTERN as2t smunu_entr[4];   //these are the sigma matrices entries
EXTERN int smunu_pos[4][6];  //and positions

//the real amd imaginary unit
EXTERN complex ONE;
EXTERN complex I;

//The base of the 16 gamma matrixes and the two rotators
EXTERN dirac_matr base_gamma[16];
EXTERN dirac_matr Pplus,Pminus;
#ifdef APPRETTO
const char gtag[16][3]={"S0","V1","V2","V3","V0","P5","A1","A2","A3","A0","T1","T2","T3","B1","B2","B3"};
#else
extern const char gtag[16][3];
#endif

