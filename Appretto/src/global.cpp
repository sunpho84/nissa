#pragma once

int **glb_coord,glb_size[4],glb_vol;
int **loc_coord,loc_size[4],loc_vol;
int *glb_of_loc_ind;
int  nproc_dir[4]={0,0,0,0};
int proc_coord[4]={0,0,0,0};
int rank,rank_tot,cart_rank;

int big_endian;

MPI_Comm cart_comm;


typedef double complex[2];

typedef complex spin[4];
typedef complex color[3];

typedef spin colorspin[3];
typedef color spincolor[4];

typedef spin spinspin[4];
typedef spinspin colorspinspin[3];
